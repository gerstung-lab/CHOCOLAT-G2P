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


class ModelGenotyping(PyroModule):
    """
    A Pyro probabilistic model class for genotyping cancer regions based on spatial transcriptomics data,
    modelling perturbation probabilities from barcode counts.

    Parameters
    ----------
    plasmid_matrix : torch.Tensor
        A matrix linking plasmids to barcodes used in the experiment.
    region_ids : torch.Tensor
        Tensor containing region identifiers for each spot in the dataset.
    data_shape : tuple of int
        Shape of the data matrix (visium spots, barcodes).
    max_range : int, optional
        The maximum number of copies of a plasmid per clone, to account for dosage-dependent variation. 
        Default is 6.

    Attributes
    ----------
    device : torch.device
        The computation device (CUDA or CPU) determined based on availability.
    n_regions : int
        Number of unique regions based on `region_ids`.
    cn_range : torch.Tensor
        Range tensor from 0 to `max_range` - 1.
    noise : torch.Tensor
        A tensor representing a fixed probability of no perturbation presence in normal tissue regions.
    data_shape : tuple of int
        Shape of the data (spots, barcodes), defining the dimensions for modeling.

    Methods
    -------
    forward(data, total_counts, mask=None)
        Defines the generative process for the observed data given the model parameters, estimating
        the spot sensitivity, barcode expression rates, and region-specific perturbation counts.
        Requires data tensor dim=(spots, barcodes) and total UMI counts per spot dim=(spots, 1)
        
    Notes
    -----
    Bayesian modeling of perturbation probabilities from barcode counts:
    The observed expression count matrix D_{s,b} (spots s by barcode genes b) is assumed to follow
    a Negative Binomial distribution with mean λ_{s,b} and overdispersion ϕ_{b}. The overdispersion 
    parameter ϕ is sampled from a Gamma distribution (shape=1000, rate=0.03), skewed towards higher 
    values to encourage the likelihood to approximate a Poisson distribution in the absence of 
    overdispersion evidence. The mean expression for each spot, λ_{s,b}, incorporates sensitivity of
    each spot, mappings of spots to clonal nodules, expected number of integrated copies of plasmid,
    linkage of plasmids to their corresponding barcodes, and barcode-specific expression rates and 
    noise:
    
    $$\\lambda_{s,g} = \\mu_{s} \\sum_{r} A_{s,r} \\sum_{g}G_{r,g} B_{g,b} k_{b} + \\xi_{b}$$

    The per-nodule plasmid integration number, modeled as an expected count of integration events, 
    describes the probability distribution over the number of plasmid copies integrated per clone,
    assuming up to `max_range` copies.
    """
    def __init__(self, plasmid_matrix, region_ids, data_shape, max_range=6):
        super().__init__()
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.plasmid_matrix = plasmid_matrix
        self.region_ids = region_ids
        self.n_regions = len(torch.unique(region_ids))
        self.cn_range = torch.tensor(np.arange(max_range), device=self.device)
        self.max_range = max_range
        self.noise = torch.ones((1, self.plasmid_matrix.shape[0], self.max_range), device=self.device) * 0.001
        self.noise[:,:,0] = 1 - 0.001 * (max_range - 1)
        self.data_shape = data_shape

    def forward(self, data, total_counts, mask=None):
        r_0 = pyro.sample('r_0', dist.Exponential(torch.tensor([1.], device=self.device)).expand([self.data_shape[1]]).to_event(1))
        epsilon = pyro.sample('epsilon', dist.Exponential(torch.tensor([1.], device=self.device)).expand([self.data_shape[1]]).to_event(1))
        
        # r_ad = 0.0001
        r_ad = pyro.sample('r_ad', dist.Exponential(torch.tensor([1.], device=self.device)).expand([self.data_shape[1]]).to_event(1))
        
        r = torch.stack([r_0*0,  epsilon], dim=0)

        g = torch.gather(r, 0, self.plasmid_matrix)

        lam = pyro.sample('spot_sensitivity', dist.Gamma(torch.tensor([3.], device=self.device),
                                                     torch.tensor([0.3], device=self.device)).expand([self.data_shape[0],1]).to_event(2))
        
        b = pyro.sample('b', dist.Gamma(torch.tensor([10.], device=self.device), torch.tensor([0.01], device=self.device)))
        
        if mask:
            total_obs = pyro.sample('total_counts',
                                    dist.Poisson(lam[mask,:] * b).to_event(2),
                                    obs=total_counts[mask,:])
        else:
            total_obs = pyro.sample('total_counts',
                                    dist.Poisson(lam * b).to_event(2),
                                    obs=total_counts)

        
        # Initialize region parameters
        region_params = pyro.sample('region_params', dist.Dirichlet(torch.ones(self.max_range, device=self.device)).expand([self.n_regions-1, self.plasmid_matrix.shape[0]]).to_event(2))
        region_params = torch.concatenate([region_params, self.noise], axis=0)
        region_params_expected = (region_params * self.cn_range.reshape(1,1,self.max_range)).sum(-1)
        # For each sample, use its associated region id to select the correct Beta parameters
        selected_params = region_params_expected[self.region_ids]

        mu = lam * (selected_params @ g) + r_ad[None,:]
        phi = pyro.sample('phi', dist.Gamma(torch.tensor([1000.], device=self.device),
                                        torch.tensor([0.03], device=self.device)).expand([1, self.data_shape[1]]).to_event(1))
        if mask:
            obs = pyro.sample("obs", dist.GammaPoisson(concentration=phi[:,:],
                                       rate=phi[:,:] / mu[mask,:]).to_event(2), obs=data[mask,:])
        else:
            obs = pyro.sample("obs", dist.GammaPoisson(concentration=phi,
                                       rate=phi / mu).to_event(2), obs=data)
            

class GuideGenotyping(PyroModule):
    """
    A Pyro variational guide for the `ModelGenotyping` class, specifying the variational distributions
    for the model's parameters. This guide is designed to approximate the posterior distribution 
    of the genotyping model using parameterized variational distributions.

    Parameters
    ----------
    plasmid_matrix : torch.Tensor
        A matrix linking plasmids to barcodes used in the experiment.
    region_ids : torch.Tensor
        Tensor containing region identifiers for each spot in the dataset.
    data_shape : tuple of int
        Shape of the data matrix (spots, barcodes).
    max_range : int, optional
        The maximum number of copies of a plasmid per clone, used for constructing Dirichlet 
        priors for region-specific plasmid copy number. Default is 6.

    Attributes
    ----------
    device : torch.device
        The computation device (CUDA or CPU) determined based on availability.
    n_regions : int
        Number of unique regions based on `region_ids`.
    cn_range : torch.Tensor
        Range tensor from 0 to `max_range` - 1.
    data_shape : tuple of int
        Shape of the data (spots, barcodes), defining the dimensions for variational parameters.

    Methods
    -------
    forward(data, total_counts, mask=None)
        Executes the variational guide, sampling from variational distributions and registering
        variational parameters with Pyro's parameter store.
        
    Notes
    -----
    This guide provides variational families for parameters like spot sensitivity (λ), overdispersion (ϕ),
    and region-specific plasmid copy numbers. It uses LogNormal distributions for positive constrained 
    parameters like spot sensitivity, and Dirichlet distributions for multinomial proportions in region
    parameters. 
    """
    def __init__(self, plasmid_matrix, region_ids, data_shape, max_range=6):
        super().__init__()
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.plasmid_matrix = plasmid_matrix
        self.region_ids = region_ids
        self.n_regions = len(torch.unique(region_ids))
        self.cn_range = torch.tensor(np.arange(max_range), device=self.device)
        self.max_range = max_range
        self.data_shape = data_shape

    def forward(self, data, total_counts, mask=None):
        # Define variational parameters
        r_loc = pyro.param('r_loc', torch.ones([self.data_shape[1], 3], device=self.device), constraint=constraints.interval(-1000,100))
        r_scale = pyro.param('r_scale',  torch.ones([self.data_shape[1], 3], device=self.device)*1e-1, constraint=constraints.positive)
        b_loc = pyro.param('b_loc', torch.tensor([7.], device=self.device), constraint=constraints.interval(1,100))
        b_scale = pyro.param('b_scale',torch.tensor([0.1], device=self.device), constraint=constraints.positive)
        b = pyro.sample('b', dist.LogNormal(b_loc, b_scale))

        # Sample r[:, 0] and epsilon separately, and ensure epsilon is positive
        epsilon = pyro.sample('epsilon', dist.LogNormal(r_loc[:, 0], r_scale[:, 0]).to_event(1))
        r_0 = pyro.sample('r_0', dist.LogNormal(r_loc[:, 1], r_scale[:, 1]).to_event(1))
        r_ad = pyro.sample('r_ad', dist.LogNormal(r_loc[:, 2], r_scale[:, 2]).to_event(1))

        # variational parameters for lamda
        lam_loc = pyro.param('lam_loc', torch.ones(self.data_shape[0],1, device=self.device), constraints.interval(-10,10))
        lam_scale = pyro.param('lam_scale', torch.ones(self.data_shape[0],1, device=self.device), constraint=constraints.positive)

        lam = pyro.sample('spot_sensitivity', dist.LogNormal(lam_loc, lam_scale).to_event(2))

        # variational parameters for phi
        phi_loc = pyro.param('phi_loc', torch.ones([1, self.data_shape[1]], device=self.device) * 100)
        phi_scale = pyro.param('phi_scale', torch.ones([1, self.data_shape[1]], device=self.device), constraint=constraints.positive)

        phi = pyro.sample('phi', dist.Gamma(phi_loc, phi_scale).to_event(1))
        alpha = pyro.param("alpha", torch.ones([self.n_regions-1, self.plasmid_matrix.shape[0], self.max_range], device=self.device),
                       constraint=constraints.positive)
        region_params = pyro.sample('region_params', dist.Dirichlet(alpha).to_event(2))
        
        
def prepare_data4genotype_models(adata_dict, reporters, plasmid_matrix,
                          region_id_column='histo_annotation_num', 
                          keys=None,
                          device='cuda'):
    """
    Prepares and transforms data from AnnData objects for use with genotyping models, 
    transferring the necessary data to the specified computational device.

    Parameters
    ----------
    adata_dict : dict of AnnData
        A dictionary of AnnData objects, each containing gene expression data and metadata 
        for CHOCOLAT-G2P 10x Visium expreiments.
    reporters : list of str
        List of gene names that serve as reporters in the plasmid_matrix.
    plasmid_matrix : numpy.ndarray
        A matrix linking plasmids to barcodes or reporters used in the experiment.
    region_id_column : str, optional
        The column name in `adata.obs` that contains region identifiers. 
        Default is 'histo_annotation_num'.
    keys : list of str, optional
        A list of keys specifying which samples in `adata_dict` to process. If None, 
        all samples in `adata_dict` will be processed.
    device : str, optional
        The device to which tensors will be transferred ('cuda' or 'cpu'). If 'cuda' is 
        specified and available, tensors will be moved to GPU. Default is 'cuda'.

    Returns
    -------
    tuple
        Returns a tuple containing four elements:
        - data: A dictionary mapping each key to a tensor of expression data for reporter genes.
        - umi: A dictionary mapping each key to a tensor of UMI counts (total counts per spot).
        - region_ids_tensor: A dictionary mapping each key to a tensor of encoded region IDs.
        - plasmid_matrix_tensor: A tensor representing the plasmid_matrix transferred to the specified device.
    """
    
    if (device =='cuda') and torch.cuda.is_available():
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')
    if keys is None:
        keys = list(adata_dict.keys())
                
    plasmid_matrix_tensor = torch.tensor(plasmid_matrix, device=device)
    print(f'Prepared data for the following samples: \n{", ".join(keys)}')
    
    data = {k: torch.tensor(v[:,reporters].X.todense(), device=device) for k, v in adata_dict.items()}
    umi = {k: torch.tensor(v.X.sum(1), device=device) for k, v in adata_dict.items()}

    # Spot to region ids
    region_ids = {k: v.obs[region_id_column].astype('object').fillna('NaN').astype('category').cat.codes for k, v in adata_dict.items()}
    region_ids_tensor = {k: torch.tensor(v.values.astype(int), device=device) for k, v in region_ids.items()}
    
    return data, umi, region_ids_tensor, plasmid_matrix_tensor

        
def train_genotype_models(data, umi, region_ids_tensor, plasmid_matrix_tensor,
                          max_range=6,
                          num_iters=10000,
                          lr=0.01,
                          keys=None, 
                          save_filename=None, 
                          return_sites=None,
                          num_samples=1000,
                          seed=42,
                          device='cuda'):
    """
    Trains genotype models for each sample specified with variational inference.
    
    Parameters
    ----------
    data : dict
        A dictionary mapping sample keys to their respective expression data tensors.
    umi : dict
        A dictionary mapping sample keys to their respective UMI count tensors.
    region_ids_tensor : dict
        A dictionary mapping sample keys to tensors of encoded region IDs.
    plasmid_matrix_tensor : torch.Tensor
        A tensor representing the plasmid matrix, used in both the model and guide.
    max_range : int, optional
        The maximum number of plasmid copy numbers per region. Default is 6.
    num_iters : int, optional
        Number of iterations to run the stochastic variational inference. Default is 10000.
    lr : float, optional
        Learning rate for the Adam optimizer. Default is 0.01.
    keys : list of str, optional
        List of sample keys to process. If None, processes all keys in `data`.
    save_filename : str, optional
        Filename to save the posterior samples as a numpy file. If None, posteriors are not saved.
    return_sites : tuple of str, optional
        Specifies the model sites from which to return posterior samples. Default includes all main sites.
    num_samples : int, optional
        Number of posterior samples to generate after training. Default is 1000.
    seed : int, optional
        Random seed for reproducibility. Default is 42.
    device : str, optional
        Specifies whether to use 'cuda' or 'cpu'. If 'cuda' is available, it will be used by default.

    Returns
    -------
    tuple
        Returns a tuple containing dictionaries for models, guides, optimizers, SVI objects, and samples:
        - models: The trained Pyro models.
        - guides: The corresponding variational guides.
        - optimizers: The Adam optimizers used during training.
        - svi_dict: The SVI objects configured for inference.
        - samples_dict: Posterior samples from the trained models (this dict is saved, if filename is specified).
    """
    if (device =='cuda') and torch.cuda.is_available():
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')
    
    if keys is None:
        keys = list(data.keys())
        
    if return_sites is None:
        return_sites=("region_params", "r_0", "epsilon",
                      "spot_sensitivity", "r_ad", "obs", "phi")
        
    # Initialize an empty dictionary to store model, guide, optimizer, and svi for each sample
    models = {}
    guides = {}
    optimizers = {}
    svi_dict = {}
    samples_dict = {} 
    
    # Loop over all adata objects and train a separate model for each
    for sample_name in keys:
        # Initialize model and guide
        models[sample_name] = ModelGenotyping(plasmid_matrix_tensor, region_ids_tensor[sample_name], data[sample_name].shape, max_range=max_range)
        guides[sample_name] = GuideGenotyping(plasmid_matrix_tensor, region_ids_tensor[sample_name], data[sample_name].shape, max_range=max_range)

        # Set up the optimizer and the inference algorithm
        optimizers[sample_name] = pyro.optim.Adam({"lr": lr})
        svi_dict[sample_name] = pyro.infer.SVI(models[sample_name], guides[sample_name], optimizers[sample_name], 
                                               loss=pyro.infer.Trace_ELBO(num_particles=3))

    # Run inference separately for each sample
    loss_lists = {sample_name: [] for sample_name in keys}
    num_iterations = num_iters
    
    torch.manual_seed(seed)
    pyro.set_rng_seed(seed)
    
    if save_filename is None:
        print('Starting model training, note that the results will not be saved. If you want to save the results, specify `save_filename`')
    else: 
        print(f'Starting model training, the restults will be saved to {save_filename}')
        
    for sample_name in keys:
        pyro.clear_param_store()
        with tqdm(total=num_iterations, desc=f"Processing {sample_name}") as pbar:
            for i in range(num_iterations):
                loss = svi_dict[sample_name].step(data[sample_name], umi[sample_name])
                loss_lists[sample_name].append(loss)
                if (i + 1) % 100 == 0:
                    pbar.set_postfix(loss=loss)
                    pbar.update(100)
            # After inference, generate samples from the trained model
            predictive = Predictive(models[sample_name], guide=guides[sample_name], num_samples=1000, return_sites=return_sites)
            samples = predictive(None, None)
            samples_dict[sample_name] = {k: v.detach().cpu().numpy() for k, v in samples.items()}
        
        

    if save_filename is not None:
        np.save(save_filename, samples_dict)
                
    return models, guides, optimizers, svi_dict, samples_dict