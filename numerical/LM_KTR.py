import numpy as np
import glob
from scipy import interpolate, optimize, integrate
import warnings
import multiprocessing as mp
from functools import partial

warnings.filterwarnings('ignore')


# Definition of functions

# Value of the loglikelihood as a function of gamma
# find mean_t directly from solving the dL/dmean_t=0 equation
def calculate_log_l(gamma, event, t_final, spline):

    # Run cumulative_hazard in parallel, change the 4 for something else if you want more cores
    p =  mp.Pool(4)
    func = partial(calculate_cumulative_hazard, gamma, spline)
    cumulative_hazard = np.array(p.map(func, t_final))
    p.close()
    
    log_hazard = calculate_log_hazard(gamma, t_final, spline)
    
    mean_t = cumulative_hazard.sum() / event.sum()
    log_l = -event.sum() * np.log(mean_t) + log_hazard[event].sum() - (1 / mean_t) * cumulative_hazard.sum()

    return -log_l

def calculate_cumulative_hazard(gamma, spline, t):

    int_Veff = integrate.quad(lambda x: np.exp(gamma * spline(x)), 0, t)[0]
    return int_Veff

def calculate_log_hazard(gamma, t_final, spline):

    Veff = spline(t_final)
    return gamma * Veff

# Actual script

## LOAD DATA

t_final=np.loadtxt('JUMPTIMES') # JUMPTIME OR FINAL SIMULATION TIME 
t_total=10000000.0           # TOTAL SIMULATION TIME
T=np.loadtxt('VMB')[:,0]            # TIME FOR ALL TRAJECTORIES (CONCATENATED)
V=np.loadtxt('VMB')[:,1]             # VMB FOR ALL TRAJECTORIES (CONCATENATED)

unique_T = np.unique(T)
unique_V = np.empty_like(unique_T)
for n_t, t in enumerate(unique_T):
    unique_V[n_t] = V[np.where(T == t)].mean()

## SPLINE FIT VMB(t)
spline = interpolate.UnivariateSpline(unique_T, unique_V, s=0)  # Update smoothing factor?
event = t_final < t_total

## LIKELIHOOD OPTIMIZATION
opt = optimize.minimize_scalar(calculate_log_l, bounds=(0.0, 1), method='bounded', args=(event, t_final, spline))
gamma = opt.x

p =  mp.Pool(4)
func = partial(calculate_cumulative_hazard, gamma, spline)
cumulative_hazard = np.array(p.map(func, t_final))
p.close()

mean_t = cumulative_hazard.sum() / event.sum()
result = np.array([1/mean_t, gamma])

## EXTRACTED GAMMA AND K
print("gamma", result[1], "k", result[0])
