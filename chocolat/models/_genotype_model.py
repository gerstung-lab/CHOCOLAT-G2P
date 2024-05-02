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
                          region_id_column='histo_annotation_num', 
                          keys=None, 
                          save_filename=None, 
                          return_sites=None,
                          num_samples=1000,
                          seed=42,
                          device='cuda'):
    
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