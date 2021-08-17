## KRAMERS TIME-DEPENDENT RATES (KTR) METHOD
### ANALYTICAL FIT OF THE SURVIVAL FUNCTION

FILES:
-- `AVE_VMB`: Maximum bias averaged over multiple simulations
-- `fit_VMB.py`: Python script to extract a and b paramters from fit of `VMB=a*log(1+b*t)`
-- `JUMPTIMES`: Simulation Jump transition times 
-- `cal_KTR.py`: KTR method using the analtyical log of VMB. Extracts gamma, rate k, and p-value.

### FIRST CALCULATE THE FIT OF THE LOGARITHMIC FUNCTION TO THE AVERAGE VMB

```bash
python3 fit_VMB.py
```

### USE FITTED a AND b PARAMETERS AS INPUT, AND INCLUDE TOTAL NUMBER OF SIMULATIONS n

```bash
python3 cal_KTR.py <n> <a> <b>
```

#### Example output:

```
gamma 0.4867346223619609 k 1.7030556148664812e-11 pval 0.4329427782593286
```

