import numpy as np

def compute_quantiles(data):
    """
    Compute quantiles for each element in a dataset.

    Parameters
    ----------
    data : array_like
        The input data for which quantiles are to be computed. This can be a list, tuple,
        or an array of numbers.

    Returns
    -------
    ndarray
        An array of quantiles corresponding to each element in the input `data`. Each element's
        quantile reflects its position in the sorted order of the original data, scaled by the
        total number of elements.

    Notes
    -----
    The quantile for a given data point is calculated as (position in sorted array + 1) / total number of data points.

    Example
    -------
    >>> data = [12, 3, 5, 11, 9]
    >>> compute_quantiles(data)
    array([1. , 0.2, 0.4, 0.8, 0.6], dtype=float32)
    """

    data = np.array(data, dtype=np.float32)
    
    # Sort the data and compute the total number of points
    sorted_data = np.sort(data)
    n = len(data)
    
    quantiles = np.zeros(n, dtype=np.float32)
    
    # Compute the quantile for each data point
    for i, value in enumerate(data):
        pos = np.where(sorted_data == value)[0][0]
        quantile = (pos + 1) / n
        
        quantiles[i] = quantile
    
    return quantiles
