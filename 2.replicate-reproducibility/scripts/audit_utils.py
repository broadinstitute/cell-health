"""
Predicting Cell Health
Gregory Way, 2019

Functions to help with analyzing and interpreting audit analyses

Import only
"""

from scipy.stats import sem, t


def get_confidence(df, col='correlation'):
    """
    Get the 95% confidence intervals of a given distribution of scores.

    Used in a pandas apply on a groupby object.
    """
    n = df.shape[1]
    m = df.loc[:, col].mean()
    std_err = sem(df.loc[:, col])
    h = std_err * t.ppf((1 + 0.95) / 2, n - 1)

    ci_low = m - h
    ci_high = m + h

    df = df.assign(correlation_ci_low=ci_low,
                   correlation_ci_high=ci_high)

    return df
