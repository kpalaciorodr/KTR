# Transition rates, survival probabilities, and quality of bias from time-dependent biased simulations

Scripts with examples for the estimation of the intrinsic rate (k_o) and the quality of bias (gamma) from time-dependent biased simulations by fitting analytically or numerically the cumulative distribution function.

## Numerical fit
In the `numerical` folder, we provide a python script for averaging and numerically fitting the maximum bias, and the survival probablity (see README in that folder). The example data is for the 2D system using y-coordinate `(k_unbias=5.6e-8)`.

## Analytical fit: 
In the `analytical` folder, we provide a python script for averaging and fitting the maximum bias using a logaritmic function. A python script is provided to fit the cumulative distribution function and perform a KS-test using the optimal fitted parameters (see README in that folder). The example data is for the CDK2-ligand unbinding using a bias-deposition time of `1ps`, and `(k_exp = 0.26 +/- 0.05 s^-1)`.

