# verify.R

##  no damping no driving force
Problem size/sensitivity dimensions: `parameter × state × time`:
```R
[1]  3  2 93
```
### Default Values:
```R
        k 
0.3678794 
v y 
0 1 
```
Analytical sensitivity of the 4th to 8th time point:
```R
[1] -0.0001862786 -0.0003448274 -0.0005518036 -0.0008071998 -0.0011110069
```
Estimated sensitivity of the 4th to 8th time point:
```R
[1] -0.0001862773 -0.0003448238 -0.0005517962 -0.0008071864 -0.0011109852
```
The error is calculated as the sum of absolute differences, normalized by the number of
time-points.
The trajectory error from `gsl odeiv2` is 9.92001e-05.
The approximate sensitivity error is 0.052472.

Loading required package: deSolve

### deSolve vs Sensitivity shifted Trajectory:
 Sum of absolute differences per time-point: 0.0363155
##  with damping and driving force
Problem size/sensitivity dimensions: `parameter × state × time`:
```R
[1]  3  2 85
```
### Default Values:
```R
        k 
0.3678794 
v y 
0 1 
```
Analytical sensitivity of the 4th to 8th time point:
```R
[1] -1.698610e-05 -3.164746e-05 -5.081105e-05 -7.446240e-05 -1.025871e-04
```
Estimated sensitivity of the 4th to 8th time point:
```R
[1] -1.700068e-05 -3.168674e-05 -5.089305e-05 -7.460989e-05 -1.028276e-04
```
The error is calculated as the sum of absolute differences, normalized by the number of
time-points.
The trajectory error from `gsl odeiv2` is 0.0147438.
The approximate sensitivity error is 0.204259.
### deSolve vs Sensitivity shifted Trajectory:
 Sum of absolute differences per time-point: 0.00127821
##  with damping but no driving force
Problem size/sensitivity dimensions: `parameter × state × time`:
```R
[1]  3  2 86
```
### Default Values:
```R
        k 
0.3678794 
v y 
0 1 
```
Analytical sensitivity of the 4th to 8th time point:
```R
[1] -1.576905e-05 -2.975076e-05 -4.811173e-05 -7.083807e-05 -9.791591e-05
```
Estimated sensitivity of the 4th to 8th time point:
```R
[1] -1.578205e-05 -2.978651e-05 -4.818721e-05 -7.097484e-05 -9.814009e-05
```
The error is calculated as the sum of absolute differences, normalized by the number of
time-points.
The trajectory error from `gsl odeiv2` is 1.29091e-05.
The approximate sensitivity error is 0.176031.
### deSolve vs Sensitivity shifted Trajectory:
 Sum of absolute differences per time-point: 0.00119699
##  with little damping but no driving force
Problem size/sensitivity dimensions: `parameter × state × time`:
```R
[1]   3   2 165
```
### Default Values:
```R
        k 
0.3678794 
v y 
0 1 
```
Analytical sensitivity of the 4th to 8th time point:
```R
[1] -5.022189e-05 -9.603336e-05 -1.565249e-04 -2.316706e-04 -3.214441e-04
```
Estimated sensitivity of the 4th to 8th time point:
```R
[1] -5.024384e-05 -9.609515e-05 -1.566570e-04 -2.319119e-04 -3.218418e-04
```
The error is calculated as the sum of absolute differences, normalized by the number of
time-points.
The trajectory error from `gsl odeiv2` is 7.37131e-05.
The approximate sensitivity error is 0.235884.

### deSolve vs Sensitivity shifted Trajectory:
Sum of absolute differences per time-point: 0.00621798


