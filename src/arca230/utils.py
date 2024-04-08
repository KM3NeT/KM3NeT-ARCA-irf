import numpy as np
from scipy.stats import chi2, poisson


def create_histogram(x_low, x_high, y):
    """
    Converts the rows of a dataframe with the bin low and high edges
    To a numpy array that can be plotted as a histogram

    Parameters:
    - x_low: dataframe column with low edges bins
    - x_high: dataframe column with high edges bins
    - y: bin values

    Returns:
    - Dataframe with column representing the number of events within the cone
    """
    x_hist = np.array([x_low, x_high]).T.flatten()
    y_hist = np.array([y, y]).T.flatten()

    return x_hist, y_hist


def limit(k, mu0, cl=0.90):
    """
    Given the background is mu0, compute the limit with corresponding confidence level
    on signal resulting from observing k events.

    Parameters:
    - k: observed events
    - mu0: number of background events
    - cl: confidence level

    Returns:
    - Limit on the signal events
    """
    q = chi2.ppf(cl, 2 * k + 2)
    return 0.5 * q - mu0


def mean_limit(mu0, confidence_level=0.90):
    """
    Given the background is mu0, compute the mean limit with corresponding confidence level

    Parameters:
    - mu0: number of background events
    - cl: confidence level

    Returns:
    - Mean limit on the signal events
    """
    c = 0

    kmax = int(max(20, mu0 + 5 * np.sqrt(mu0)))

    for k in range(kmax):
        l = limit(k, mu0, confidence_level)
        c += poisson.pmf(k, mu0) * l

    return c


def nobs_disc(mu0, alpha, beta):
    """
    Given the background is mu0, compute the least detected number of signal events (mu_lds)
    for a significance and power

    Parameters:
    - mu0: number of background events
    - alpha: significance p-value
    - beta: power

    Returns:
    - Minimum number of signal events for alpha discovery
    """
    nmax = 100000
    n_crit = None

    for n_obs in range(nmax):
        pvalue = 1 - poisson.cdf(n_obs, mu=mu0)
        if pvalue < alpha:
            n_crit = 1 + n_obs
            break

    for cont in range(nmax):
        mu_lds = cont * 0.1  # optimize step size
        p = 1 - poisson.cdf(n_crit - 1, mu=mu0 + mu_lds)
        if p > (1 - beta):
            break

    return mu_lds
