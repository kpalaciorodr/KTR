import numpy as np
import glob
from scipy import interpolate, optimize, integrate
import warnings
warnings.filterwarnings('ignore')


# Definition of functions

# Value of the loglikelihood as a function of gamma
# find mean_t directly from solving dL/dmean_t=0
def calculate_log_l(gamma):
    global event, t_final, spline
    cum_hazard = calculate_cum_hazard(gamma)
    log_hazard = calculate_log_hazard(gamma)
    mean_t = cum_hazard.sum() / event.sum()
    log_l = -event.sum() * np.log(mean_t) + log_hazard[event].sum() - (1 / mean_t) * cum_hazard.sum()
    return -log_l

def calculate_cum_hazard(gamma):
    global t_final, spline
    int_Veff = np.zeros(t_final.shape[0])
    for n, t in enumerate(t_final):
        int_Veff[n] = integrate.quad(lambda x: np.exp(gamma * spline(x)), 0, t)[0]
    return int_Veff

def calculate_log_hazard(gamma):
    global t_final, spline
    Veff = np.zeros(t_final.shape[0])
    Veff = spline(t_final)
    return gamma * Veff

# Actual script

## LOAD DATA

t_final=np.load('JUMPTIMES') # JUMPTIME OR FINAL SIMULATION TIME 
t_total=10000000.0           # TOTAL SIMULATION TIME
T=np.load('TIME')            # TIME FOR ALL TRAJECTORIES (CONCATENATED)
V=np.load('VMB')             # VMB FOR ALL TRAJECTORIES (CONCATANATED)
unique_T = np.unique(T)
unique_V = np.empty_like(unique_T)

for n_t, t in enumerate(unique_T):
    unique_V[n_t] = V[np.where(T == t)].mean()

## SPLINE FIT VMB(t)
spline = interpolate.UnivariateSpline(unique_T, unique_V, s=0)  # Update smoothing factor?
event = t_final < t_total

## LIKELIHOOD OPTIMIZATION
opt = optimize.minimize_scalar(calculate_log_l, bounds=(0.0, 1), method='bounded')
gamma = opt.x
mean_t = calculate_cum_hazard(gamma).sum() / event.sum()
result = np.array([1/mean_t, gamma])

## EXTRACTED GAMMA AND K
print("gamma", result[1],"k", result[0])
