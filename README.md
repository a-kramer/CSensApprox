# CSensApprox
Experiments with Sensitivity Approximation (SA) in C using the gsl ode solvers.
This program will solve an ODE model for parameters stored in an hdf5 file and any input vector stored in the same hdf5 file.

The depenencies should be minimal, for now this is:

1. GNU Scientific Library
2. HDF5 and specifically HDF5LT

The model will be loaded from a shared library, to keep the model subject to change.

## Manual

Example usage:

```bash
MODEL="SomeModel" # with ./${MODEL}.so in `pwd`
./gsl_odeiv -m ${MODEL} -d ${MODEL}.h5 -t "-90:1:60" --abs-tol 1e-8 --rel-tol 1e-6 1> ${MODEL}_out.tsv
```

Options:

|short|`--`long|value type|default|meaning|
|:----:|:---|:--------:|:-----:|:------|
|`-a` |`--abs-tol`| double | 10⁻⁶ | absolute tolerance (solver)|
|`-r` |`--rel-tol`| double | 10⁻⁵ | relative tolerance (solver)|
|`-m` |`--model`| string | NONE | name of the model (no suffix)|
|`-d` |`--data` | string | NONE | simulation instructions (hdf5 file)|
|`-h` |`--step-size`| double | 10⁻³ | initial step size (suggestion) |
|`-t` |`--sim-time`| double | NONE | solver [time span](#Simulation-Time-Span) | 

The time span is a range (e.g `-1:0.1:6`), with an increment (explained in another section).

### Simulation Instructions (HDF5)

The `-d` option provides the solver with model setups. The file needs
a `prior` group with a `mu` dataset with values for the parameters of
the model (µ is in log-space). Additionally, there needs to be a
`data` group, which contains datasets (with any name) with some data
values not used here. The data should come with hdf5 ATTRIBUTES:

| attribute name | meaning |
|---------------:|:--------|
| time | measurement times for the data|
| input | input parameters in linear space (`u`) |
| InitialValue | the initial value for the state vector `y` (for each simulation run) |

For this program, the values in the hdf5 DATASET are not important,
only the simulation instructions are used. The maximum of the _time_
vector is used to set a reasonable simulation time-span.

The model has a parameter slot (a vector `p`), which will be assembled like this:

```octave
p=[exp(mu), u]; # 
```

This is because the model has unknown parameters, for which we have a probability density with median µ (logarithmic) and known parameters u, which represent the input to the model (as part of an _experiment_).

### Tolerances and Step-Size

The tolerances inform the step size controls, and are used as described in the [gsl documentation](https://www.gnu.org/software/gsl/doc/html/ode-initval.html#adaptive-step-size-control):

```octave
Dᵢ = abs_tol + rel_tol × |yᵢ|
```

where `Dᵢ` is the integration error bound. The `step-size` is adaptive, this is only a suggestion.


### Simulation Time-Span

The time span can be given in GNU Octave/Matlab/Julia syntax: `a:b:c`
or alternatively, just as `'a b c'`. This will be interpreted as 

a. initial value
b. increment
c. final value

The parser is very simple, so `a:c` will not work as expected. If
anything is 0, the value will be corrected using some defaults and
assumptions. The default for `c` is `1.3 * max(measurement time)` in
the hdf5 data file.
