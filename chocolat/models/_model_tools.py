import numpy as np

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
