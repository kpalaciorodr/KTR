## KRAMERS TIME-DEPENDENT RATES (KTR) METHOD
### NUMERICAL FIT OF THE SURVIVAL FUNCTION

FILES:
* `VMB`: Simulation time `t` and Maximum bias `VMB(t)` for multiple simulations (concatenated).
         Text format: (first column) time [steps], (second column) maximum bias [kT].
* `FINALTIMES`: Final simulation time (for non-events) or time to jump (for events)
               for multiple simulations.
* `LM_KTR.py`: Python script to numerically fit the averaged maximum bias, survival probability and optimize
               the likelihood function.

### DEPENDENCIES

Please have installed in the python environment: `scipy`, `pandas`, and `numpy`.

### CALCULATE THE FIT OF THE MAXIMUM BIAS, SURVIVAL FUNCTION AND OPTIMIZE THE LIKELIHOOD.

The total simulation time t_total is an input argument:

```bash
python3 LM_KTR.py <t_total>
```

#### Example output:

```
> python3 LM_KTR.py 10000000
gamma 0.8016868722270777 k 5.5145856981263585e-08
KS TEST Statistics 0.06937221515955472 Prob 0.6952255696129455
```
