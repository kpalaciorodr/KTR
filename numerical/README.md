## KRAMERS TIME-DEPENDENT RATES (KTR) METHOD
### NUMERICAL FIT OF THE SURVIVAL FUNCTION

FILES:
* `VMB`: Maximum bias `VMB(t)` for multiple simulations. Numpy format (concatenated)
* `TIME`: Simulation time `t` for multiple simulations. Numpy format (concatenated)
* `JUMPTIMES`: Final simulation time (for non-events) or time to jump (for events)
              for multiple simulation. (Numpy format)
* `LM_KTR.py`: Python script to fit numericaly the survival function and optimize
              the likelihood function
              
### DEPENDENCIES

### CALCULATE THE FIT OF THE SURVIVAL FUNCTION AND OPTIMIZE THE LIKELIHOOD

```bash
python3 LM_KTR.py
```

#### Example output:

```
gamma 0.7863379248552048 k 5.620335920122176e-08
```

