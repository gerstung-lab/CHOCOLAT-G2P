import pyro
from pyro.infer.autoguide import AutoDiagonalNormal
from pyro.nn import PyroSample
from pyro.distributions import constraints
from pyro.nn import PyroModule
from pyro.infer import SVI, TraceMeanField_ELBO, Predictive
from pyro.distributions.transforms.ordered import OrderedTransform
import pyro.distributions as dist
import torch
import numpy as np
from tqdm import tqdm
import pandas as pd
import os

from ._model_tools import get_sample_prior_parametrisation

class ModelGlmmBernoulli(PyroModule):
    """
    A Generalised Linear Mixed Model (GLMM) for binary outcomes with multiple predictors, implemented in Pyro. Example of a binary outcome is IHC status of the nodule region.
    This model infers the probability of expression for each gene (IHC) across different plasmids and samples,
    accounting for region-specific perturbation probabilities and sample-specific biases.
    Note that data is supplied on a per region basis.

    Parameters
    ----------
    inferred_spot_sensitivity : numpy.ndarray
        Inferred sensitivity of each spot (visium spots) from genotype model.
    inferred_region_params : numpy.ndarray
        Inferred region parameters (plasmid, regions) from genotype model.
    region_ids : numpy.ndarray
        Array containing region identifiers for each spot in the dataset.
    sample_indices : numpy.ndarray
        Array containing sample identifiers for each spot in the dataset.
    data_shape : tuple of int
        Shape of the expression data matrix (genes, plasmids).

    Attributes
    ----------
    device : torch.device
        The computation device (CUDA or CPU) determined based on availability.
    inferred_spot_sensitivity : torch.Tensor
        Tensor of inferred spot sensitivities.
    inferred_region_params : torch.Tensor
        Tensor of inferred region parameters. 
    data_shape : torch.Tensor
        Tensor of the shape of the data (genes, plasmids).
    region_ids : torch.Tensor
        Tensor of region identifiers.
    sample_indices : torch.Tensor
        Tensor of sample identifiers.
    n_samples : int
        Number of unique samples in the dataset.
    n_genes : int
        Number of genes in the dataset.
    n_plasmids : int
        Number of plasmids in the dataset.
    n_regions : int
        Number of unique regions in the inferred region parameters.

    Methods
    -------
    forward(expression_data)
        Defines the generative process for the binary expression data given the model parameters,
        estimating the probability of expression for each gene across plasmids and samples.
    """

    def __init__(self, inferred_spot_sensitivity, inferred_region_params, region_ids, sample_indices, data_shape):
        super().__init__()
        device = torch.device("cpu" if torch.cuda.is_available() else "cpu")
        self.device = device
        print(f'using {device}') 
        self.inferred_spot_sensitivity = torch.tensor(inferred_spot_sensitivity, device=device)
        self.inferred_region_params = torch.tensor(inferred_region_params, device=device)
        self.data_shape = torch.tensor(np.array(data_shape), device=device)
        self.region_ids = torch.tensor(np.array(region_ids).astype(int), device=device)
        self.sample_indices = torch.tensor(np.array(sample_indices).astype(int), device=device)
        self.n_samples = np.max(sample_indices) + 1
        self.n_genes = data_shape[0]
        self.n_plasmids = data_shape[1]
        self.n_regions = self.inferred_region_params.shape[1]
        
    def forward(self, expression_data):

        region_params = pyro.sample('region_params', 
                                    dist.Beta(self.inferred_region_params[:,:,0], self.inferred_region_params[:,:,1]).to_event(2))

        weights = pyro.sample('w', dist.Laplace(torch.tensor([0.], device=self.device),
                                                torch.tensor([0.1], device=self.device)).expand((self.data_shape[0],
                                                                        self.data_shape[1] + 1)).to_event(2))
        

        sample_bias = pyro.sample('sb', dist.Laplace(torch.tensor([0.], device=self.device),
                                                torch.tensor([1.], device=self.device)).expand((1, self.n_samples)).to_event(2))
        
        
        
        linear_term = (weights[:,:self.data_shape[1]]  @ (1-region_params[:,:]) + weights[:,self.data_shape[1]:]) + sample_bias[:,self.sample_indices] 
        
        lam = 1./ (1. + torch.exp(- linear_term))
        pyro.sample('obs', dist.Bernoulli(lam).to_event(2), obs=expression_data)


class GuideGlmmBernoulli(PyroModule):
    """
    A Pyro variational guide for the `GLMM_Binary_mult` class, specifying the variational distributions
    for the model's parameters. This guide is designed to approximate the posterior distribution 
    of the GLMM model using parameterised variational distributions. Note that data is supplied on a per region basis.

    Parameters
    ----------
    inferred_spot_sensitivity : numpy.ndarray
        Inferred sensitivity of each spot (visium spots) from genotype model.
    inferred_region_params : numpy.ndarray
        Inferred region parameters (plasmid, regions) from genotype model.
    region_ids : numpy.ndarray
        Array containing region identifiers for each spot in the dataset.
    sample_indices : numpy.ndarray
        Array containing sample identifiers for each spot in the dataset.
    data_shape : tuple of int
        Shape of the expression data matrix (genes, plasmids).

    Attributes
    ----------
    device : torch.device
        The computation device (CUDA or CPU) determined based on availability.
    data_shape : torch.Tensor
        Tensor of the shape of the data (genes, plasmids).
    n_samples : int
        Number of unique samples in the dataset.
    inferred_spot_sensitivity : torch.Tensor
        Tensor of inferred spot sensitivities.
    inferred_region_params : torch.Tensor
        Tensor of inferred region parameters.
    n_spots : int
        Number of spots (visium spots) in the dataset.
    n_genes : int
        Number of genes in the dataset.
    n_plasmids : int
        Number of plasmids in the dataset.
    n_regions : int
        Number of unique regions in the inferred region parameters.
    total_params : int
        Total number of variational parameters, including sample biases and weights.

    Methods
    -------
    forward(expression_data)
        Executes the variational guide, sampling from variational distributions and registering
        variational parameters with Pyro's parameter store.
    """
    def __init__(self, inferred_spot_sensitivity, inferred_region_params, region_ids, sample_indices, data_shape):
        super().__init__()
        device = torch.device("cpu" if torch.cuda.is_available() else "cpu")
        self.device = device
        self.data_shape = torch.tensor(np.array(data_shape), device=device)
        self.n_samples = np.max(sample_indices) + 1
        self.inferred_spot_sensitivity = torch.tensor(inferred_spot_sensitivity, device=device)
        self.inferred_region_params = torch.tensor(inferred_region_params, device=device)
        self.n_spots = inferred_spot_sensitivity.shape[0]
        self.n_samples = np.max(sample_indices) + 1
        self.n_genes = data_shape[0]
        self.n_plasmids = data_shape[1]
        self.n_regions = self.inferred_region_params.shape[1]
        
        self.total_params = self.data_shape[0] * (self.data_shape[1] + 1) + self.n_samples


    def forward(self, expression_data):
        
        mean_param = pyro.param('mean_param', torch.zeros(self.total_params, device=self.device))
        sigma_param = pyro.param('sigma_scale', torch.ones(self.total_params, device=self.device) * 0.05,
                                 constraint=constraints.positive)

        joint_dist = dist.Normal(mean_param, sigma_param)
        
        samples = pyro.sample('joint', joint_dist.to_event(1), infer={'is_auxiliary': True})
        sb_sample = samples[:self.n_samples].reshape(1, -1)
        w_sample = samples[self.n_samples:].reshape(self.data_shape[0], self.data_shape[1]+1)
        pyro.sample('sb', dist.Delta(sb_sample).to_event(2))
        pyro.sample('w', dist.Delta(w_sample).to_event(2))
        spot_loc = pyro.param('spot_loc', torch.zeros(self.n_spots, device=self.device), constraint=constraints.interval(-5, 5))
        spot_scale = pyro.param('spot_scale', torch.ones(self.n_spots, device=self.device) * 0.05, constraint=constraints.positive)
        
class ModelGlmmPoisson(PyroModule):
    """
    A Generalised Linear Mixed Model (GLMM) for count data using a Poisson distribution.
    This model infers the expected count of expression for each gene across different plasmids and samples,
    accounting for region-specific perturbation probabilities and sample-specific biases. Note that data is supplied on a per spot basis.

    Parameters
    ----------
    inferred_spot_sensitivity : numpy.ndarray
        Inferred sensitivity of each spot (visium spots) from a previous model.
    inferred_region_params : numpy.ndarray
        Inferred region parameters (plasmid, regions) from a previous model.
    region_ids : numpy.ndarray
        Array containing region identifiers for each spot in the dataset.
    sample_indices : numpy.ndarray
        Array containing sample identifiers for each spot in the dataset.
    data_shape : tuple of int
        Shape of the expression data matrix (genes, visium spots).

    Attributes
    ----------
    device : torch.device
        The computation device (CUDA or CPU) determined based on availability.
    inferred_spot_sensitivity : torch.Tensor
        Tensor of inferred spot sensitivities.
    inferred_region_params : torch.Tensor
        Tensor of inferred region parameters.
    data_shape : torch.Tensor
        Tensor of the shape of the data (genes, plasmids).
    region_ids : torch.Tensor
        Tensor of region identifiers.
    sample_indices : torch.Tensor
        Tensor of sample identifiers.
    n_samples : int
        Number of unique samples in the dataset.

    Methods
    -------
    forward(expression_data)
        Defines the generative process for the count expression data given the model parameters,
        estimating the expected count of expression for each gene across plasmids and samples.
    """

    def __init__(self, inferred_spot_sensitivity, inferred_region_params, region_ids, sample_indices, data_shape):
        super().__init__()
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.device = device
        print(f'using {device}') 
        self.inferred_spot_sensitivity = torch.tensor(inferred_spot_sensitivity, device=device)
        self.inferred_region_params = torch.tensor(inferred_region_params, device=device)
        self.data_shape = torch.tensor(np.array(data_shape), device=device)
        self.region_ids = torch.tensor(np.array(region_ids).astype(int), device=device)
        self.sample_indices = torch.tensor(np.array(sample_indices).astype(int), device=device)
        self.n_samples = np.max(sample_indices) + 1
        
    def forward(self, expression_data):
        spot_sensitivity = pyro.sample('spot_sensitivity', 
                                       dist.LogNormal(self.inferred_spot_sensitivity[:,0], self.inferred_spot_sensitivity[:,1]).to_event(1))

        region_params = pyro.sample('region_params', 
                                    dist.Beta(self.inferred_region_params[:,:,0], self.inferred_region_params[:,:,1]).to_event(2))

        weights = pyro.sample('w', dist.Laplace(torch.tensor([0.], device=self.device),
                                                torch.tensor([0.001]*8 + [1], device=self.device)).expand((self.data_shape[0],
                                                                                                self.data_shape[1]+1)).to_event(2))
        sample_bias = pyro.sample('sb', dist.Laplace(torch.tensor([0.], device=self.device),
                                                torch.tensor([1.], device=self.device)).expand((1, self.n_samples)).to_event(2))
        linear_term = (weights[:,:self.data_shape[1]]  @ (1-region_params[:,self.region_ids]) + weights[:,self.data_shape[1]:]) + sample_bias[:,self.sample_indices] 
        clipped_linear_term = torch.clamp(linear_term, min=-8, max=8)
        lam = spot_sensitivity[None,:] * torch.exp(clipped_linear_term)
        pyro.sample('obs', dist.Poisson(lam).to_event(2), obs=expression_data)


class GuideGlmmPoisson(PyroModule):
    """
    A Pyro variational guide for the `ModelGlmmPoisson` class, specifying the variational distributions
    for the model's parameters. This guide is designed to approximate the posterior distribution 
    of the GLMM model using parameterised variational distributions.

    Parameters
    ----------
    inferred_spot_sensitivity : numpy.ndarray
        Inferred sensitivity of each spot (visium spots) from a previous model.
    inferred_region_params : numpy.ndarray
        Inferred region parameters (plasmid, regions) from a previous model.
    region_ids : numpy.ndarray
        Array containing region identifiers for each spot in the dataset.
    sample_indices : numpy.ndarray
        Array containing sample identifiers for each spot in the dataset.
    data_shape : tuple of int
        Shape of the expression data matrix (genes, plasmids).

    Attributes
    ----------
    device : torch.device
        The computation device (CUDA or CPU) determined based on availability.
    data_shape : torch.Tensor
        Tensor of the shape of the data (genes, plasmids).
    n_samples : int
        Number of unique samples in the dataset.
    inferred_spot_sensitivity : torch.Tensor
        Tensor of inferred spot sensitivities.
    inferred_region_params : torch.Tensor
        Tensor of inferred region parameters.
    n_spots : int
        Number of spots (visium spots) in the dataset.
    total_params : int
        Total number of variational parameters, including sample biases and weights.

    Methods
    -------
    forward(expression_data)
        Executes the variational guide, sampling from variational distributions and registering
        variational parameters with Pyro's parameter store.
    """
    
    def __init__(self, inferred_spot_sensitivity, inferred_region_params, region_ids, sample_indices, data_shape):
        super().__init__()
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.device = device
        self.data_shape = torch.tensor(np.array(data_shape), device=device)
        self.n_samples = np.max(sample_indices) + 1
        self.inferred_spot_sensitivity = torch.tensor(inferred_spot_sensitivity, device=device)
        self.inferred_region_params = torch.tensor(inferred_region_params, device=device)
        self.n_spots = inferred_spot_sensitivity.shape[0]
        self.total_params = self.data_shape[0] * (self.data_shape[1] + 1) + self.n_samples

    def forward(self, expression_data):
        mean_param = pyro.param('mean_param', torch.zeros(self.total_params, device=self.device))
        sigma_param = pyro.param('sigma_scale', torch.ones(self.total_params, device=self.device) * 0.05,
                                 constraint=constraints.positive)
        joint_dist = dist.Normal(mean_param, sigma_param)
        samples = pyro.sample('joint', joint_dist.to_event(1), infer={'is_auxiliary': True})
        sb_sample = samples[:self.n_samples].reshape(1, -1)
        w_sample = samples[self.n_samples:].reshape(self.data_shape[0], self.data_shape[1]+1)
        pyro.sample('sb', dist.Delta(sb_sample).to_event(2))
        pyro.sample('w', dist.Delta(w_sample).to_event(2))
        spot_loc = pyro.param('spot_loc', torch.zeros(self.n_spots, device=self.device), constraint=constraints.interval(-5, 5))
        spot_scale = pyro.param('spot_scale', torch.ones(self.n_spots, device=self.device) * 0.05, constraint=constraints.positive)
        spot_sensitivity = pyro.sample('spot_sensitivity', 
                       dist.LogNormal(spot_loc, spot_scale).to_event(1))

def load_ihc_annotation(adata_dict, binary_feature_list, keys,
                        folder='../data/IHCs_12ROIs/'):
    """
    Loads immunohistochemistry (IHC) annotation data and integrates it into the AnnData objects.

    Parameters
    ----------
    adata_dict : dict
        Dictionary of AnnData objects, indexed by sample keys.
    binary_feature_list : list of str
        List of binary feature names to be loaded and processed.
    keys : list of str
        List of sample keys corresponding to the data files.
    folder : str, optional
        Path to the folder containing the IHC annotation files. Default is '../data/IHCs_12ROIs/'.

    Returns
    -------
    observations_dict : dict
        Dictionary containing the processed binary feature annotations for each sample.
    """

    observations_dict = {}

    for key in keys:
        temp_dict = {} 
        for val_keys in binary_feature_list:
            adata_dict[key].obs[f'{val_keys}_stain'] = pd.read_csv(os.path.join(folder,
                                                                                f'{key}_{val_keys}.csv'),
                                                            index_col=0).loc[adata_dict[key].obs.index,:]
            adata_dict[key] = _spread_ihc_annotations(adata_dict[key],
                                                     f'{val_keys}_stain')
            temp_dict[val_keys] = adata_dict[key].obs[f'spread_{val_keys}_stain'].notna().astype(int)

        observations_dict[key] = pd.DataFrame(temp_dict)
    return observations_dict
        
        
def prepare_data4glmm(adata_dict, samples_dict, observations_dict,
                             region_id_column='histo_annotation_num', 
                             keys=None):
    """
    Prepares data for Generalised Linear Mixed Models (GLMM) by reindexing and stacking region and sample identifiers.

    Parameters
    ----------
    adata_dict : dict
        Dictionary of AnnData objects, indexed by sample keys.
    samples_dict : dict
        Dictionary containing sample information.
    observations_dict : dict
        Dictionary containing binary feature annotations for each sample.
    region_id_column : str, optional
        Column name in AnnData objects representing region identifiers. Default is 'histo_annotation_num'.
    keys : list of str, optional
        List of sample keys to be processed. If None, all keys in `adata_dict` are used.

    Returns
    -------
    prior_parametrisation_tensor_dict : dict
        Dictionary of tensors containing prior parametrisation information for each sample.
    region_indices_stacked : numpy.ndarray
        Array of stacked region indices.
    sample_indices_by_region : numpy.ndarray
        Array of sample indices averaged by region.
    sample_indices_stacked : numpy.ndarray
        Array of stacked sample indices.
    observations_tensor : torch.Tensor
        Tensor containing the processed observations data.
    observations_by_region : torch.Tensor
        Tensor of observations averaged by region.
    """

    if keys is None:
        keys = list(adata_dict.keys())
        
    prior_parametrisation_dict = {} 
           
    for key in keys:
        prior_parametrisation_dict[key] = get_sample_prior_parametrisation(samples_dict[key])
        
    region_ids = {k: v.obs[region_id_column].astype('object').fillna('NaN').astype('category').cat.codes.copy() for k, v in adata_dict.items()}
    region_ids_isnormal = {k: (v.obs[region_id_column].astype('object').fillna('NaN').astype('category') == 'NaN').copy() for k, v in adata_dict.items()}

    for k in region_ids.keys():
        region_ids[k][region_ids_isnormal[k]] = -1    
    
    new_region_ids, back_mapping, sample_mapping = _reindex_region_ids(region_ids, keys)
    
    (sample_indices_stacked,
     region_indices_stacked,
     observations_tensor,
     prior_parametrisation_tensor_dict) = _stack_data_tensors(new_region_ids,
                                                             observations_dict,
                                                             prior_parametrisation_dict,
                                                             keys)

    filter_ids = np.where(np.array(region_indices_stacked) != -1)[0]

    sample_indices_stacked = np.array(sample_indices_stacked)[filter_ids]
    region_indices_stacked = np.array(region_indices_stacked)[filter_ids]
    observations_tensor = observations_tensor[filter_ids]
    prior_parametrisation_tensor_dict['spot_sensitivity'] = prior_parametrisation_tensor_dict['spot_sensitivity'][filter_ids,:]
    
    observations_by_region = np.zeros((region_indices_stacked.max() + 1, observations_tensor.shape[1]))
    sample_indices_by_region = torch.zeros(region_indices_stacked.max() + 1)
    for i in range(observations_by_region.shape[0]):
        observations_by_region[i,:] = np.array(observations_tensor.cpu().numpy()[region_indices_stacked == i, :]).mean(0).astype(int)
        sample_indices_by_region[i] = int(sample_indices_stacked[region_indices_stacked == i].mean())
    
    sample_indices_by_region = np.array(sample_indices_by_region).astype(int)
    observations_by_region = torch.tensor(observations_by_region)

    return prior_parametrisation_tensor_dict, region_indices_stacked, sample_indices_by_region, sample_indices_stacked, observations_tensor, observations_by_region

def train_glm_models(prior_parametrisation_tensor_dict,
                     region_indices,
                     sample_indices,
                     observations_tensor,
                     obs,
                     model,
                     guide,
                     num_iters=1000,
                     lr=0.01,
                     seed=42,
                     device='cuda'):
    """
    Trains Generalised Linear Mixed Models (GLMM) using Stochastic Variational Inference (SVI).

    Parameters
    ----------
    prior_parametrisation_tensor_dict : dict
        Dictionary of tensors containing prior parametrisation information for each sample.
    region_indices : numpy.ndarray
        Array of region indices.
    sample_indices : numpy.ndarray
        Array of sample indices.
    observations_tensor : torch.Tensor
        Tensor containing the processed observations data (visium spots, genes).
    obs : torch.Tensor
        Tensor of observed data (the acutal data that model is trained on).
    model : PyroModule
        GLMM model class to be trained.
    guide : PyroModule
        Variational guide class for the model.
    num_iters : int, optional
        Number of iterations for the SVI optimisation. Default is 1000.
    lr : float, optional
        Learning rate for the Adam optimiser. Default is 0.01.
    seed : int, optional
        Random seed for reproducibility. Default is 42.
    device : str, optional
        Device for computation ('cuda' or 'cpu'). Default is 'cuda'.

    Returns
    -------
    mean_marginal : numpy.ndarray
        Mean of the marginal posterior distributions of the model weights.
    variance_marginal : numpy.ndarray
        Variance of the marginal posterior distributions of the model weights.
    loss_list : list of float
        List of loss values for each iteration of the SVI optimisation.
    """

    if (device =='cuda') and torch.cuda.is_available():
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')
                    
    model = model(prior_parametrisation_tensor_dict['spot_sensitivity'],
            prior_parametrisation_tensor_dict['region_params'],
            region_indices,
            sample_indices,
           [observations_tensor.shape[1], 8])

    guide = guide(prior_parametrisation_tensor_dict['spot_sensitivity'],
                prior_parametrisation_tensor_dict['region_params'],
                region_indices,
                sample_indices,
               [observations_tensor.shape[1], 8])

    optimizer = pyro.optim.Adam({"lr": lr})
    svi = pyro.infer.SVI(model, guide,  optimizer, loss=pyro.infer.Trace_ELBO(num_particles=3))

    torch.manual_seed(seed)
    pyro.set_rng_seed(seed)

    num_iterations = num_iters
    
    loss_list = []

    # move to device
    obs_tensor = obs.to(model.device)
#     observations_by_region = observations_by_region.to(model.device)

    pyro.clear_param_store()
    with tqdm(total=num_iterations, desc=f"Processing") as pbar:
        for i in range(num_iters):
            loss = svi.step(obs_tensor.T)
            loss_list.append(loss)
            if (i + 1) % 100 == 0:
                pbar.set_postfix(loss=loss) 
                pbar.update(100) 
                
    mean_param = pyro.param('mean_param').detach().cpu().numpy().copy()
    sigma_param = pyro.param('sigma_scale').detach().cpu().numpy().copy()


    mean_marginal = mean_param[(model.n_samples):].reshape(model.data_shape[0]
                                                , model.data_shape[1]+1)

    variance_marginal = sigma_param[(model.n_samples):].reshape(model.data_shape[0]
                                                , model.data_shape[1]+1)

                
    return mean_marginal, variance_marginal, loss_list

def _reindex_region_ids(region_ids, keys):
    """
    Reindexes region IDs to ensure unique identifiers across samples.

    Parameters
    ----------
    region_ids : dict
        Dictionary of region IDs, indexed by sample keys.
    keys : list of str
        List of sample keys corresponding to the data files.

    Returns
    -------
    new_region_ids : dict
        Dictionary of reindexed region IDs.
    back_mapping : dict
        Dictionary mapping new region IDs to original region IDs.
    sample_mapping : dict
        Dictionary mapping new region IDs to sample keys.
    """

    new_region_ids = {}
    back_mapping = {}
    sample_mapping = {}
    max_region_id = -1  # Initialise with -1 to get the true maximum in the first iteration
    
    for sample in keys:
        regions = region_ids[sample]
        # Shift the region IDs by max_region_id + 1, but keep -1 as is
        shifted_regions = regions.apply(lambda x: x + max_region_id + 1 if x != -1 else x)
        
        # Update the maximum region ID based on the shifted values
        max_region_id = shifted_regions.max()
        
        # Store the reindexed region IDs
        new_region_ids[sample] = shifted_regions
        
        # Create back-mapping and sample-mapping for each shifted region ID
        for old, new in zip(regions, shifted_regions):
            if old != -1:  # Do not include -1, as it is common across samples
                back_mapping[new] = old
                sample_mapping[new] = sample
                
    return new_region_ids, back_mapping, sample_mapping


def _stack_data_tensors(region_ids, observations_dict, prior_parametrisation_dict, keys):
    """
    Stacks data tensors for regions, observations, and prior parametrisations.

    Parameters
    ----------
    region_ids : dict
        Dictionary of region IDs, indexed by sample keys.
    observations_dict : dict
        Dictionary containing binary feature annotations for each sample.
    prior_parametrisation_dict : dict
        Dictionary of prior parametrisation information for each sample.
    keys : list of str
        List of sample keys corresponding to the data files.

    Returns
    -------
    tuple
        Contains stacked sample indices, region indices, observations tensor, and prior parametrisation tensors.
    """

    sample_indices = []
    region_indices = []
    observations = []
    prior_parametrisation_tensor_dict = {'region_params': [], 'spot_sensitivity': []}
    
    for sample_index, sample_key in enumerate(keys):
        regions = region_ids[sample_key]
        observation = observations_dict[sample_key]
        region_param = prior_parametrisation_dict[sample_key]['region_params']
        spot_sensitivity = prior_parametrisation_dict[sample_key]['spot_sensitivity']
        
        n_observations = len(regions)
        
        sample_indices.extend([sample_index] * n_observations)
        region_indices.extend(regions)
        
        observations.append(observation)
        prior_parametrisation_tensor_dict['region_params'].append(np.array(region_param).T)
        prior_parametrisation_tensor_dict['spot_sensitivity'].append(np.array(spot_sensitivity).T)
    observations_tensor = torch.tensor(np.concatenate(observations, axis=0).astype(np.float32))
    prior_parametrisation_tensor_dict['region_params'] = torch.tensor(np.concatenate(prior_parametrisation_tensor_dict['region_params'], axis=1).astype(np.float32))
    prior_parametrisation_tensor_dict['spot_sensitivity'] = torch.tensor(np.concatenate(prior_parametrisation_tensor_dict['spot_sensitivity']).astype(np.float32))
    return sample_indices, region_indices, observations_tensor, prior_parametrisation_tensor_dict

def _spread_ihc_annotations(adata, key_name):
    """
    Spreads IHC annotations by assigning the most abundant annotation that exceeds a 15% threshold
    to spots within each histopathological annotation.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    key_name : str, optional
        Name of the IHC annotation key. For example 'IHC-GFP'.

    Returns
    -------
    AnnData
        Annotated data matrix with spread IHC annotations.
    """

    adata_obs = adata.obs.copy()
    filtered_data = adata_obs.dropna(subset=['histo_annotation', key_name])

    grouped_counts = filtered_data.groupby(['histo_annotation', key_name]).size().reset_index(name='counts')

    total_counts = grouped_counts.groupby('histo_annotation')['counts'].transform('sum')
    grouped_counts['proportion'] = grouped_counts['counts'] / total_counts * 100

    filtered_grouped_counts = grouped_counts[grouped_counts['proportion'] > 15]
    most_abundant_annotation = filtered_grouped_counts.groupby('histo_annotation').apply(lambda x: x.nlargest(1, 'counts')).reset_index(drop=True)

    annotation_to_gfp = {row['histo_annotation']: row[key_name] for _, row in most_abundant_annotation.iterrows()}

    adata_obs['spread_' + key_name] = adata_obs['histo_annotation'].map(annotation_to_gfp).fillna(adata_obs[key_name])
    adata.obs = adata_obs
    return adata


