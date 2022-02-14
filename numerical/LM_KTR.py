#####NUMERICAL FIT OF VMB
## RUN:
## python3 LM_KTR.py <t_total>

### INPUTS AND PARAMETERS
# t_total = total simulation time
# t = vector with final simulation time (for non-events) or jump time (for events)


import numpy as np
from sys import argv
from scipy import interpolate, optimize, integrate, stats
import warnings
import multiprocessing as mp
from functools import partial
import pandas as pd

warnings.filterwarnings("ignore")


# Definition of functions

# Value of the loglikelihood as a function of gamma
# find mean_t directly from solving the dL/dmean_t=0 equation
def calculate_log_l(gamma, event, t, spline):

    # Run cum_hazard in parallel, change the 4 for something else if you want more cores
    p = mp.Pool(4)
    func = partial(calculate_cum_hazard, gamma, spline)
    cum_hazard = np.array(p.map(func, t))
    p.close()

    log_hazard = calculate_log_hazard(gamma, t, spline)

    mean_t = cum_hazard.sum() / event.sum()
    log_l = -event.sum() * np.log(mean_t) + log_hazard[event].sum() - (1 / mean_t) * cum_hazard.sum()

    return -log_l


def calculate_cum_hazard(gamma, spline, t):

    int_Veff = integrate.quad(lambda x: np.exp(gamma * spline(x)), 0, t)[0]
    return int_Veff


def calculate_log_hazard(gamma, t, spline):

    Veff = spline(t)
    return gamma * Veff


# Compute CDF for our model using fitted parameters
def cumulative_density_fct(gamma, mean_t, spline, t):
    # Run cum_hazard in parallel, change the 4 for something else if you want more cores
    p = mp.Pool(4)
    func = partial(calculate_cum_hazard, gamma, spline)
    cum_hazard = np.array(p.map(func, t))
    p.close()
    return 1 - np.exp(-cum_hazard / mean_t)


# Perform one sample  Kolmogorov-Smirnov test
def ks_1samp_test(gamma, mean_t, spline, event, t):
    # Only take times with events
    N = np.sum(event)
    x = np.sort(t)[:N]
    cdfvals = cumulative_density_fct(gamma, mean_t, spline, x)
    Dplus = (np.arange(1.0, N + 1) / N - cdfvals).max()
    Dminus = (cdfvals - np.arange(0.0, N) / N).max()
    D = np.max([Dplus, Dminus])
    prob = stats.distributions.kstwo.sf(D, N)
    prob = np.clip(prob, 0, 1)
    return D, prob  # {"statistic": D, "pvalue": prob}


def main():
    # Load data and parameters
    t = np.loadtxt("FINALTIMES")
    t_total = int(argv[1])  # TOTAL SIMULATION TIME
    T = np.loadtxt("VMB")[:, 0]  # TIME FOR ALL TRAJECTORIES (CONCATENATED)
    V = np.loadtxt("VMB")[:, 1]  # VMB FOR ALL TRAJECTORIES (CONCATENATED)

    # Calculating average VMB
    unique_T = np.unique(T)
    df = pd.DataFrame({"T": T, "V": V})
    unique_V = df.groupby("T").mean()["V"].values

    # SPLINE FIT VMB(t)
    spline = interpolate.UnivariateSpline(unique_T, unique_V, s=0,ext=1)
    event = t < t_total

    # LIKELIHOOD OPTIMIZATION
    opt = optimize.minimize_scalar(calculate_log_l, bounds=(0.0, 1), method="bounded", args=(event, t, spline))
    gamma = opt.x

    p = mp.Pool(4)
    func = partial(calculate_cum_hazard, gamma, spline)
    cum_hazard = np.array(p.map(func, t))
    p.close()

    mean_t = cum_hazard.sum() / event.sum()
    result = np.array([1 / mean_t, gamma])
    ks_stat, ks_prob = ks_1samp_test(gamma, mean_t, spline, event, t)
    # EXTRACTED GAMMA AND K AS WELL AS KS-TEST RESULT
    print("gamma", result[1], "k", result[0])
    print("KS TEST Statistics", ks_stat, "Prob", ks_prob)


if __name__ == "__main__":
    main()
