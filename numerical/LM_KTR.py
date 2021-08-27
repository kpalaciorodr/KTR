import numpy as np
import glob
from scipy import interpolate, optimize, integrate
import warnings
import multiprocessing as mp
from functools import partial

warnings.filterwarnings('ignore')


# Definition of functions

def calculate_log_l(gamma, event, t_final, spline):

    # Run cum_hazard in parallel, change the 4 for something else if you want more cores
    p =  mp.Pool(4)
    func = partial(calculate_cum_hazard, gamma, spline)
    cum_hazard = np.array(p.map(func, t_final))
    p.close()
    
    log_hazard = calculate_log_hazard(gamma, t_final, spline)
    
    mean_t = cum_hazard.sum() / event.sum()
    log_l = -event.sum() * np.log(mean_t) + log_hazard[event].sum() - (1 / mean_t) * cum_hazard.sum()

    return -log_l

def calculate_cum_hazard(gamma, spline, t):

    int_Veff = integrate.quad(lambda x: np.exp(gamma * spline(x)), 0, t)[0]
    return int_Veff

def group(key, value):
    """
    group the values by key
    returns the unique keys, and the average of the values per-key
    By: Eelco Hoogendoorn
    https://stackoverflow.com/questions/7790611/average-duplicate-values-from-two-paired-lists-in-python-using-numpy
    """

    # Upcast to numpy arrays
    key = np.asarray(key)
    value = np.asarray(value)

    # First, sort by key
    I = np.argsort(key)
    key = key[I]
    value = value[I] # Sort values according to how the keys were sorted

    # The slicing points of the bins to sum over, i.e, where the sorted keys stop being unique
    slices = np.concatenate(([0], np.where(key[:-1] != key[1:])[0]+1))

    # First entry of each bin is a unique key
    unique_keys = key[slices]

    # Sum over the slices specified by index
    per_key_sum = np.add.reduceat(value, slices)

    # Number of counts per key is the difference of our slice points. cap off with number of keys for last bin
    key_count = np.diff(np.append(slices, len(key)))

    # Calculate the mean for the values
    mean_values = per_key_sum/key_count
    
    return unique_keys, mean_values

def calculate_log_hazard(gamma, t_final, spline):

    Veff = spline(t_final)
    return gamma * Veff

# Actual script

## LOAD DATA

t_final=np.load('JUMPTIMES') # JUMPTIME OR FINAL SIMULATION TIME 
t_total=10000000.0           # TOTAL SIMULATION TIME
T=np.load('TIME')            # TIME FOR ALL TRAJECTORIES (CONCATENATED)
V=np.load('VMB')             # VMB FOR ALL TRAJECTORIES (CONCATENATED)

unique_T, unique_V = group(T, V)

## SPLINE FIT VMB(t)
spline = interpolate.UnivariateSpline(unique_T, unique_V, s=0)  # Update smoothing factor?
event = t_final < t_total

## LIKELIHOOD OPTIMIZATION
opt = optimize.minimize_scalar(calculate_log_l, bounds=(0.0, 1), method='bounded', args=(event, t_final, spline))
gamma = opt.x

p =  mp.Pool(4)
func = partial(calculate_cum_hazard, gamma, spline)
cum_hazard = np.array(p.map(func, t_final))
p.close()

mean_t = cum_hazard.sum() / event.sum()
result = np.array([1/mean_t, gamma])

## EXTRACTED GAMMA AND K
print("gamma", result[1],"k", result[0])
