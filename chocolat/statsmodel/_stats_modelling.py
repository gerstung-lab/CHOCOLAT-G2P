import numpy as np
import scipy.stats as st
from tqdm import tqdm

def genotype_prob(plasmid_probs_beta_params, genotype_matrix):
    """
    Computes the likelihood of genotypes based on plasmid probabilities using Beta-binomial distribution.

    Parameters
    ----------
    plasmid_probs_beta_params : numpy.ndarray
        Array containing the alpha and beta parameters of the Beta distribution for plasmid probabilities (2, n_nodules, n_plasmids).
    genotype_matrix : numpy.ndarray
        Matrix containing genotypes, where each row represents a genotype and each column corresponds to a plasmid.

    Returns
    -------
    numpy.ndarray
        Array of likelihoods for each genotype.
    """

    prob_1 = st.betabinom(n=1,
             a=plasmid_probs_beta_params[0],
             b=plasmid_probs_beta_params[1]).pmf(0)
    likely_candidates = []
    for i in range(genotype_matrix.shape[0]):
        prob_genotype = (np.array([1,0])[genotype_matrix[i]].reshape(-1,1) + np.array([-1,1])[genotype_matrix[i]].reshape(-1,1) * prob_1.T)
        prob_prod = np.prod(prob_genotype, 0)
        likely_candidates.append(prob_prod)
    return np.array(likely_candidates)

def resample_population(plasmid_probs_beta_params, n_samples=5000):
    """
    Resamples the population of genotypes based on plasmid probabilities using Beta-binomial distribution.

    Parameters
    ----------
    plasmid_probs_beta_params : numpy.ndarray
        Array containing the alpha and beta parameters of the Beta distribution for plasmid probabilities (2, n_nodules, n_plasmids).
    n_samples : int, optional
        Number of samples to generate. Default is 5000.

    Returns
    -------
    tuple
        - numpy.ndarray: Probabilities of plasmids being present.
        - numpy.ndarray: Simulated genotypes.
    """

    prob_1 = st.betabinom(n=1,
         a=plasmid_probs_beta_params[0],
         b=plasmid_probs_beta_params[1]).pmf(0)
    sim_genotypes = np.array([np.random.binomial(1,np.repeat(prob_1[None,:,i], n_samples, axis=0)) for i in range(8)])
    return prob_1, sim_genotypes

def marginal_frequency(sim_genotypes, axis=0, rep_axis=1, bins=9, brange=None, desnity=True):
    """
    Computes the marginal frequency of simulated genotypes.

    Parameters
    ----------
    sim_genotypes : numpy.ndarray
        Simulated genotypes. E.g. (n_plasmids, n_samples, n_nodules)
    axis : int, optional
        Axis along which to compute the marginal frequency. Default is 0.
    rep_axis : int, optional
        Repetition axis. Default is 1.
    bins : int, optional
        Number of bins for histogram. Default is 9.
    brange : tuple, optional
        Range for histogram bins. Default is None.
    density : bool, optional
        Whether to normalize the histogram. Default is True.

    Returns
    -------
    list
        List of marginal frequencies.
    """

    ps = []
    np.transpose(sim_genotypes, (axis, rep_axis, 3 - (axis + rep_axis)))
    if brange is None:
        brange = (0, bins-1)
    for i in tqdm(range(sim_genotypes.shape[1])):
        p, b = np.histogram(sim_genotypes.sum(0)[i], bins=bins, range=brange, density=True)
        ps.append(p)
    return ps

def compute_odds_ratio_plasmids(prob_1):
    """
    Computes the odds ratio matrix for plasmids based on their presence probabilities.

    Parameters
    ----------
    prob_1 : numpy.ndarray
        Probabilities of plasmids being present (n_nodules, n_plasmids).

    Returns
    -------
    numpy.ndarray
        Log odds ratio matrix.
    """

    dim = prob_1.shape[-1]
    i = 0
    j = 0
    log_odds_matrix = np.zeros((dim, dim))

    prob_success = prob_1[:,::-1]

    for i in range(dim):
        for j in range(dim):
            p11 = prob_success[:, i] * prob_success[:, j]
            p01 = (1 - prob_success[:, i]) * prob_success[:, j]
            p10 = prob_success[:, i] * (1 - prob_success[:, j])
            p00 = (1 - prob_success[:, i]) * (1 - prob_success[:, j])

            log_odds_matrix[i,j] = np.log10(p11.mean() * p00.mean() / (p01.mean() * p10.mean()))

    return log_odds_matrix


def compute_observed_genotypes_probability(prob_1, ref_genotype_table, n_samples=5000, bins=25):
    """
    Computes the expected occurrence of observed genotypes based on plasmid probabilities.

    Parameters
    ----------
    prob_1 : numpy.ndarray
        Probabilities of plasmids being present (n_nodules, n_plasmids).
    ref_genotype_table : numpy.ndarray
        Matrix of reference genotype combinations (n_genotypes, n_plasmids).
    n_samples : int, optional
        Number of samples to generate. Default is 5000.
    bins : int, optional
        Number of bins for histogram. Default is 25.

    Returns
    -------
    numpy.ndarray
        Array of expected occurrences for observed genotypes.
    """

    expected_occurance_observed = []
    for i in tqdm(range(ref_genotype_table.shape[0])):
        prob_genotype = (np.array([1,0])[ref_genotype_table[i]].reshape(-1,1) + np.array([-1,1])[ref_genotype_table[i]].reshape(-1,1) * prob_1.T)
        prob_prod = np.prod(prob_genotype, 0)

        sim = np.random.binomial(1,np.repeat(prob_prod[None,:], n_samples, axis=0))
        x, bins = np.histogram(sim.sum(1),bins=bins, range=(0,bins))
        expected_occurance_observed.append(x/n_samples)
        
    return np.array(expected_occurance_observed)
    
def hdpi_sim_genotypes_counts(prob_array, n_samples=1000, p=5):
    """
    Computes the highest density posterior intervals (HDPI) for simulated genotypes counts.

    Parameters
    ----------
    prob_array : numpy.ndarray
        Array of probabilities for each genotype (n_genotypes, n_count_registers).
    n_samples : int, optional
        Number of samples to generate. Default is 1000.
    p : float, optional
        Percentile for HDPI. Default is 5.

    Returns
    -------
    numpy.ndarray
        Array of HDPI for each genotype.
    """

    hdpis = []
    for i in range(prob_array.shape[0]):
        sim = np.argmax(np.random.multinomial(1, prob_array[i], size=n_samples), axis=1)
        hdpis.append([np.median(sim), np.percentile(sim, p/2), np.percentile(sim, 100 - p/2)])
    return np.array(hdpis)

def simulate_hnull_genotype_probability(prob_1, ref_genotype_table, n_samples=5000, bins=25):
    """
    Simulates the null genotype probabilities based on observed plasmid probabilities.

    Parameters
    ----------
    prob_1 : numpy.ndarray
        Probabilities of plasmids being present (n_nodules, n_plasmids).
    ref_genotype_table : numpy.ndarray
        Reference table of genotype combinations (n_genotypes, n_plasmids).
    n_samples : int, optional
        Number of samples to generate. Default is 5000.
    bins : int, optional
        Number of bins for histogram. Default is 25.

    Returns
    -------
    numpy.ndarray
        Expected null genotype probabilities.
    """

    plasmid_probs = prob_1.mean(0)
    plasmid_probs = plasmid_probs / plasmid_probs.sum()

    scaled_probs = plasmid_probs * np.mean(prob_1.sum(1))

    n_observations = prob_1.shape[0]
    
    hnull_samples = np.random.binomial(np.repeat([np.repeat([[1]* prob_1.shape[1]], n_observations, 0)], n_samples, 0),
                                   np.repeat([np.repeat([list(scaled_probs)], n_observations, 0)], n_samples, 0))
    
    # Convert each 8-element binary vector to a decimal number
    decimals = np.array([[_binary_to_decimal(row) for row in array] for array in hnull_samples])
    # Count the frequency of each number (0-255) for each of the 1000 entries
    frequency_matrix = np.array([np.bincount(entry, minlength=2**prob_1.shape[1]) for entry in decimals])
    frequency_matrix_sorted = frequency_matrix[:,[_binary_to_decimal(x) for x in ref_genotype_table]]
    
    expected = []
    for i in range(2**prob_1.shape[1]):
        p_x, bins = np.histogram(frequency_matrix_sorted[:,i],bins=bins, range=(0,bins))
        p_x = p_x / n_samples
        expected.append(p_x)
    return np.array(expected)

def dif_sim_genotypes(observed_genotype_prob, expected_genotype_prob, n_samples=1000):
    """
    Computes the difference in simulated genotype counts between observed and expected probabilities.

    Parameters
    ----------
    observed_genotype_prob : numpy.ndarray
        Observed genotype probabilities (n_genotypes, n_count_registers).
    expected_genotype_prob : numpy.ndarray
        Expected genotype probabilities (n_genotypes, n_count_registers).
    n_samples : int, optional
        Number of samples to generate. Default is 1000.

    Returns
    -------
    list
        List of differences in simulated genotype probabilities.
    """

    ar_bigger = []
    for i in tqdm(range(observed_genotype_prob.shape[0])):
        ar_exp = []
        ar_obs = []
        for n in range(n_samples):
            ar_exp.append(np.argmax(np.random.multinomial(1, np.array(expected_genotype_prob).T[:,i])))
            ar_obs.append(np.argmax(np.random.multinomial(1, np.array(observed_genotype_prob).T[:,i])))
        ar_bigger.append((((np.array(ar_obs) > np.array(ar_exp)).sum()) + ((np.array(ar_obs) == np.array(ar_exp)).sum()) * 0.5) / n_samples)
    return ar_bigger

def _binary_to_decimal(binary_vector):
    """
    Converts a binary vector to a decimal number.

    Parameters
    ----------
    binary_vector : array-like
        Binary vector to be converted.

    Returns
    -------
    int
        Decimal representation of the binary vector.
    """

    return int("".join(str(int(bit)) for bit in binary_vector), 2)

