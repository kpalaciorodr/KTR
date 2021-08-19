## KRAMERS TIME-DEPENDENT RATES (KTR) METHOD
### ANALYTICAL FIT OF THE SURVIVAL FUNCTION

FILES:
* `AVE_VMB`: Maximum bias averaged over multiple simulations. Text format: (first column) time [ps], (second column) averaged maximum bias.
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
>python3 fit_VMB.py
a 6.752938658295208 b 0.025656565586767554
>python3 cal_KTR.py 50 6.752938658295208 0.025656565586767554
gamma 0.4869417571263078 k 1.678147638540419e-11 pval 0.5118174752056448
```

With k in ps^-1
