# The Harmonic Oscillator

This example demonstrates that the sensitivity approximation is not
tied to applications in systems biology, nor is it required to use the
SBtab format to create the necessary hdf5 files.

The simulation instructions
[HarmonicOscillator.h5](./HarmonicOscillator.h5) have been written
using the R program
[MakeHarmonicOscillatorData.R](./MakeHarmonicOscillatorData.R). The
datasets themselves are unimportant, as described in the main
[README.md](../README.md). 

The model itself is in [HarmonicOscillator.vf](HarmonicOscillator.vf),
a [vfgen](https://warrenweckesser.github.io/vfgen/) file that was used
to automatically create both the
[gsl](https://www.gnu.org/software/gsl/doc/html/index.html) source
[HamonicOscillator_gvf.c](HarmonicOscillator.c) and the R source
[HarmonicOscillator.R](HarmonicOscillator.R) for the deSolve
package. The auto generated [demo](HarmonicOscillator_demo.R) shows
how simulate the model in R.

## Initial Value Problem

The ODE, with function `y(x)` and independent variable `x` has this form:

We use `m=1` (unit mass):
```
           y''= -k*y - c*y',   y(0) = y0,    y'(0) = v0
```
where `y'` is the first derivative `dy(x)/dx` and `y''` the second derivative. This equation leads to oscillation and the solution is better parameterized using `w=sqrt(k)` (angular velocity) and `r`, the damping ratio: `r = c/2w`.

## Solution for the Autonomous System with no Driving Force

Here, we assume the there is no external driving force `F=0` (in the model file).

```
y = a*exp(-r*w*x)*cos(sqrt(1-r*r)*w*x+f)
```

Check:
```
y' = -a*exp(-r*w*x) * (r*w) * cos(sqrt(1-r*r)*w*x+f)
     -a*exp(-r*w*x)         * sin(sqrt(1-r*r)*w*x+f) * sqrt(1-r*r)*w
y' = - y*(r*w)
     - a*exp(-r*w*x)*sin(sqrt(1-r*r)*w*x+f)*sqrt(1-r*r)*w

y'' = - y'*r*w
      + a*exp(-r*w*x)*(r*w)*sin(sqrt(1-r*r)*w*x+f)*sqrt(1-r*r)*w
      - y * (1-r*r)*w*w

y'' = - y'*r*w
      + (-y*r*r*w*w - y'*r*w)
      - y * (1-r*r)*w*w

y'' = - 2*y'*r*w
      - y * r*r*w*w
      - y*ww
      + y * r*r*w*w

y'' = -y'*c - y*k      [success]
                ^
              (k=w*w)
```

How should the amplitude `a` and phase `f` be chosen to match the initial values of the state variables?

```
y(0) = + a*cos(f) = y0
v(0) = - a*cos(f)*r*w - a*sin(f)*sqrt(1-r*r)*w = v0
     = - y0*r*w - y0*tan(f)*sqrt(1-r*r)*w
     
=> -(v0 + y0*r*w)/(y0*sqrt(1-r*r)*w) = tan(f)
```
and consequently 
```
   f = atan(-(v0 + y0*r*w)/(y0*sqrt(1-r*r)*w))
   a = y0/cos(f)
```
The sensitivity `S(x;k)` of this solution is just the derivative with repsect to `k`: `dy(x;k)/dk`. In R this results in:
```R
srr1 <- sqrt(1-rr)
dwdk <- 1/(2*w)
	
S <- - a * exp(-r*w*x) * (r*dwdk*x) * cos(srr1*w*x + f) 
     - a * exp(-r*w*x) *              sin(srr1*w*x + f) * (srr1*dwdk*x)
```

# Comparisons

The R script `verify.R` makes 3 different comparisons: 
1. Compare the `gsl_odeiv2` calculated trajectory to the exact solution (always with with `F=0`)
2. Perform a trajectory shift using the sensitivity matrix: `y(x;k+dk)
   ?= y(x;k) + S(x;k)*dk`, and compare the shifted trajector to the exact
   solution at `k+dk`, where `dk` is a small additive shift.
3. Perform the same trajectory shift, but compare to a numerical
   solution at `k+dk`, using the deSolve package
   
The third comparison covers cases, for which the analytical solution
is ill-equipped (constant driving forces).

It produces some text output and figures of these comparisons. The
text output from one run follows from here.


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


