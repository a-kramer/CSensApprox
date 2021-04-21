# CaMKIIs

This model was created from the [SBtab](sbtab.net) files in this
repository: [a-kramer/CaMKIIs](https://github.com/a-kramer/CaMKIIs).

The model simlates the time courses of molecular reactants that are
involved in a network of biochemical reactions. The simulations are
grouped in different experiments, each _major experiment_ is a _dose
reponse curve_. Each curve contains several points, each representing
the real biochemical systems response to a constant stimulation (the
_dose_). These individual responses are _minor experiments_: each
minor experiment requires a simulation of the model to reproduce,
because each data point is the final point of a time course created by
applying a different input to the model.

In aggregate, all curves contain 99 points, so there are 99
simulations that need to be performed.

This example also illustrates how to interpret the results using core
[MATLAB](https://www.mathworks.com/products/matlab.html). 

It is mostly possible to use [GNU
Octave](https://www.gnu.org/software/octave/index), but it lacks the
hdf5 functions to read hdf5 attributes. Octave imports _all_ contents of
an hdf5 file in one go, omitting all attributes. The results of
`gsl_odeiv` simulations don't have any attributes, only the data file
does.

## Input and Output

We subdivide the model's parameters into two groups: the unknown
parameters _ρ_ and the known parameters _u_. The known parameters are
related to the input to the model, a different phrasing is that
together with the initial values they reproduce the _experimental
conditions_ (e.g. a treatment dose). 

Usually, many reaction rate coefficients (and similar biological
parameters) are unknown or have large confidence intervals. So, we
assume that the parameters are to be found through optimization or
Bayesian sampling. 

For this reason, each dataset in the hdf5 data file comes with an
`input`, `InitialValue` attribute, a and `time` attribute (measurement
time) as well as some indexing information (major index, minor index, and
overall _flat_ index). The indexing is to preserve the order the user
intended.

The data file also includes the parameters _µ_ and _σ_ (or _Σ_) of a
log-normal prior distribution for the unknown parameters.

The ODE model interface has a slot for `t` (time), `y` (state) and `p` (parameters of any kind). 

So, we combine the different classes of parameters into one vector: 
```matlab
p = [exp(mu),u]
```

## Verification of Results

We don't have an analytical solution for this model. Instead we test
the sensitivity, as we did in the other examples

- shift the trajectory, using the estimated sensitivity: _y+Δy = y(t;p) + S(t;p)*Δp_
- simulate the shifted trajectory directly via `ode15s`: _y(t;p+Δp)_
- take the norm of the difference betweeen the predicted trajectory _y+Δdy_ and the numerical solution _y(t;dy)_
   + calculate the norm of the difference
   + normalize by the norm of the sum
   
The matlab files [matlab_verify.m](./matlab_verify.m) and
[rel_err.m](./rel_err.m) perform this comparison.
