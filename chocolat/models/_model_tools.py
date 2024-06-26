import numpy as np
import anndata
import pandas as pd
import scanpy as sc

def _estimate_beta_params(samples):
    """
    Estimates the alpha and beta parameters of a Beta distribution using the method of moments.

    Parameters
    ----------
    samples : array-like
        Sample data for which the Beta distribution parameters are estimated.

    Returns
    -------
    alpha_est : float
        Estimated alpha parameter.
    beta_est : float
        Estimated beta parameter.
    """

    samples = np.array(samples)
    # Calculate the sample mean
    sample_mean = np.mean(samples, axis=0)
    
    # Calculate the sample variance
    sample_var = np.var(samples, axis=0)
    
    # Use the method of moments to estimate the alpha and beta parameters
    common_factor = (sample_mean * (1 - sample_mean) / sample_var) - 1
    alpha_est = sample_mean * common_factor
    beta_est = (1 - sample_mean) * common_factor
    
    return alpha_est, beta_est

def get_sample_prior_parametrisation(samples):
    """
    Computes the prior parametrisation for a given sample using Beta distribution parameters
    and the log-normal distribution for spot sensitivity.

    Parameters
    ----------
    samples : dict
        Dictionary containing sample data with keys 'region_params' and 'spot_sensitivity'.

    Returns
    -------
    dict
        Dictionary with estimated parameters for 'region_params' and 'spot_sensitivity'.
    """

    alpha, beta = _estimate_beta_params(samples['region_params'][:,:,:,0])
    
    mu = np.log(samples['spot_sensitivity'][:,:,0]).mean(0)
    sigma = np.log(samples['spot_sensitivity'][:,:,0]).std(0)
    
    return {'region_params': [alpha, beta], 'spot_sensitivity': [mu, sigma]}

def get_model_weights_array(mean_marginal, variance_marginal, gene_list, remove_genes=None, sort=True):
    """
    Processes and filters the GLM model weights array based on mean and variance, optionally removes specified genes, and sorts the data.

    Parameters
    ----------
    mean_marginal : numpy.ndarray
        Array of mean marginal posterior distributions of the model parameters (n_genes, n_plasmids+1).
    variance_marginal : numpy.ndarray
        Array of variance marginal posterior distributions of the model parameters (n_genes, n_plasmids+1).
    gene_list : list of str
        List of gene names corresponding to the model parameters.
    remove_genes : list of str, optional
        List of genes to be removed from the results. Default is None.
    sort : bool, optional
        Whether to sort the data using UMAP-based dimensionality reduction. Default is True.

    Returns
    -------
    data : numpy.ndarray
        Filtered and optionally sorted array of model weights.
    data_names : numpy.ndarray
        List of gene names corresponding to the filtered and optionally sorted model weights.
    """

    test = np.zeros_like(mean_marginal)
    test[np.sign(mean_marginal) > 0] = (mean_marginal[np.sign(mean_marginal) > 0] -
                                        3*variance_marginal[np.sign(mean_marginal) > 0]**0.5) > 0
    test[np.sign(mean_marginal) <= 0] = (mean_marginal[np.sign(mean_marginal) <= 0] +
                                        3*variance_marginal[np.sign(mean_marginal) <= 0]**0.5) < 0
    if remove_genes is not None:
        test[np.isin(gene_list, remove_genes)] = False
    test = test.astype(bool)
    mean_marginal_altered = mean_marginal.copy()

    filter_array = test[:,:-1].max(1)

    data = mean_marginal_altered[filter_array, :-1].copy()
    data_names = np.array(gene_list)[filter_array]
    
    if sort == True:
        glm_pred_ad = anndata.AnnData(X=data, obs=pd.DataFrame(data_names).set_index(0))
        sc.pp.neighbors(glm_pred_ad, n_neighbors=30, n_pcs=None)
        sc.tl.umap(glm_pred_ad, min_dist=1, spread=2,n_components=1)
        umap_1d_coords = glm_pred_ad.obsm['X_umap'].copy()

        leaf_order = np.argsort(umap_1d_coords[:,0])

        data = data[leaf_order][:,:].T[:,:]
        data_names = list(glm_pred_ad[leaf_order].obs_names)
        
    return data, data_names