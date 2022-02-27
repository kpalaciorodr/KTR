## KRAMERS TIME-DEPENDENT RATES (KTR) METHOD
### NUMERICAL FIT OF THE SURVIVAL FUNCTION

FILES:
* `VMB`: Simulation time t and Maximum bias VMB(t) for multiple simulations.
         Text format: (first column) time [ps], (second column) maximum bias [kT].
* `FINALTIMES`: Final simulation time (for non-events) or time to jump (for events)
               for multiple simulations.
* `LM_KTR.py`: Python script to numerically fit the averaged maximum bias, survival probability and optimize
               the likelihood function.

### DEPENDENCIES

Please have installed in the python environment: `scipy`, `pandas`, and `numpy`.

### CALCULATE THE FIT OF THE MAXIMUM BIAS, SURVIVAL FUNCTION AND OPTIMIZE THE LIKELIHOOD.

The total simulation time t_total is an input argument:

```bash
python3 LM_KTR.py <t_total> <VMB_file> <FINALTIMES_file>
```

#### Example output:

```
> python3 LM_KTR.py 40000000 VMB FINALTIMES
gamma 0.7633359239575043 k 6.052093359190397e-08
```
