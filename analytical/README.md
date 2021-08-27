## KRAMERS TIME-DEPENDENT RATES (KTR) METHOD
### ANALYTICAL FIT OF THE SURVIVAL FUNCTION

FILES:
* `VMB`:  Simulation time `t` and Maximum bias `VMB(t)` for multiple simulations (concatenated). 
          Text format: (first column) time [ps], (second column) maximum bias [kT].
* `fit_VMB.py`: Python script to extract a and b paramters from fit of `VMB=a*log(1+b*t)`
* `JUMPTIMES`: Simulation-jump transition times. Text format: (first column) simulation-jump time [ps].
* `cal_KTR.py`: KTR method using the analtyical log of VMB. Extracts gamma, rate k, and p-value.

### DEPENDENCIES

Please have installed in the python environment: `scipy` and `numpy`.

### FIRST CALCULATE THE FIT OF THE LOGARITHMIC FUNCTION TO THE AVERAGE VMB

```bash
python3 fit_VMB.py
```

### USE FITTED a AND b PARAMETERS AS INPUT, AND INCLUDE TOTAL NUMBER OF SIMULATIONS n

```bash
python3 cal_KTR.py <n> <a> <b>
```

#### Example output:

```bash
> python3 fit_VMB.py
a 6.625614290702432 b 0.02724910306967592
> python3 cal_KTR.py 50 6.625614290702432 0.02724910306967592
gamma 0.48895491848631023 k 1.746803196680984e-11 pval 0.5006084116663858
```

With k in ps^-1
