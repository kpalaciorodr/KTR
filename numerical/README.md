## KRAMERS TIME-DEPENDENT RATES (KTR) METHOD
### NUMERICAL FIT OF THE SURVIVAL FUNCTION

FILES:
* `VMB`: Simulation time `t` and Maximum bias `VMB(t)` for multiple simulations (concatenated)
* `JUMPTIMES`: Final simulation time (for non-events) or time to jump (for events)
               for multiple simulations. 
* `LM_KTR.py`: Python script to fit numericaly the survival function and optimize
               the likelihood function
              
### DEPENDENCIES

Please have installed in the python environment: `scipy` and `numpy`.

### CALCULATE THE FIT OF THE SURVIVAL FUNCTION AND OPTIMIZE THE LIKELIHOOD

```bash
python3 LM_KTR.py
```

#### Example output:

```
gamma 0.8016868722318178 k 5.5145856980254623e-08
```

