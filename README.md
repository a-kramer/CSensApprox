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
|`-m` |`--model`| string | first `*.so` file | name of the model (no suffix)|
|`-d` |`--data` | string | `${MODEL}.h5` | simulation instructions (hdf5 file)|
|`-h` |`--step-size`| double | 10⁻³ | initial step size (suggestion) |
|`-t` |`--sim-time`| double | 0:1%:max(tⱼ) | solver _time span_ | 

If no model name is given, the program searches the current working
directory for the first file ending in `.so` and takes everything
before `.so` as the model's name.

The [time span](#Simulation-Time-Span) is a range, e.g: `-1:0.1:6`, with an increment.

### Simulation Instructions (HDF5)

We use hdf5 files, because almost any language can read and write hdf5
files (julia, R, matlab and python have packages with bindings to the
C API that work very well). These hdf5 files are well suited for the transfer
of data between different architectures. 

But, one thing that hdf5 is not that well suited for is preserving the
indexing order of an array that has more than one index.

To hdf5 the order of the indices is from slow to fast: `(i,j,k)` means
`i*Nj*Hk+j*Nk+k`. From a certain point of view this is the same order
as an array in C would have, if it is indexed like this:
`a[i][j][k]`. But, C arrays are rarely indexed like this (because it
is not easy to pass such a thing to a function), and usually we would
do the index arithmetic ourselves: e.g. `A[i+j*Ni+k*Nj*Ni]`. So,
really the decision what the fast index is in C is up to the
individual programmer.

The R language has the reverse ordering of indices: `[i,j,k]` means
`i+j*Ni+k*Nj*Ni`.  So, don't expect the arrays to have the same order
of dimensions in R compared to what `h5dump` or `h5ls` report.

In any case the memory underneath is unaffected and it is always
possible to preserve the content by re-ordering the size information,
at no point do we need to transpose anything in memory.

The `-d` option provides the solver with model setups in hdf5
form. The file needs a `prior` group (a group is a bit like a folder)
with a `mu` dataset with values for the parameters of the model (µ is
interpreted as being in log-space). Additionally, there needs to be a
`data` group, which contains datasets (any name will do) with some
data values. We don't actually use these values here, this content is
intended for parameter estimation where the model parameters are
adjusted to fit this data. But, to make a simulation of these datasets
possible, simulation instructions must be provided. We do this in the
form of hdf5 ATTRIBUTES:

| attribute name | meaning |
|---------------:|:--------|
| time | measurement times for the data|
| input | input parameters in linear space (`u`) |
| InitialValue | the initial value for the state vector `y` (for each simulation run) |
| index | a number from 0 to n-1, given n datasets |

For this program (it just solves the initial value problem once), the
values in the DATASET are not important, only the simulation
instructions are used. The maximum of the _time_ vector is used to set
a reasonable simulation time-span.

Writing various datasets to a file will not preserve the order. So, to
keep the order of the datasets as the user intended them, we also
include an index attribute. We use it to import these instructions
exactly into that order (and simulate them in that order as well).

The model has a parameter slot (a vector `p`), which will be assembled like this:

```R
p=c(exp(mu), u); # in R, for illustration
```

This is because the model is assumed to have unknown parameters, for
which we have a probability density with median _µ_ (log-normal
distribution) and known parameters _u_, which represent the input to
the model (as part of an _experiment_). In C, this concatenation is of
course more difficult to write.

### Tolerances and Step-Size

The tolerances inform the step size controls, and are used as described in the [gsl documentation](https://www.gnu.org/software/gsl/doc/html/ode-initval.html#adaptive-step-size-control):

```C
D[i] = abs_tol + rel_tol * |y[i]|
```

where `D[i]` is the integration error bound. The `step-size` is adaptive, this is only a suggestion.


### Simulation Time-Span

The time span can be given in GNU Octave/Matlab/Julia syntax: `a:b:c`
or alternatively, just as `'a b c'`. This will be interpreted as 

|symbol|meaning        |default     |
|-----:|:--------------|-----------:|
|    a | initial value | 0.0        |
|    b | increment     | 1%         |
|    c | final value   | 1.3 max(t) |

The parser is very simple, so `a:c` will not work as expected (this
will be interpreted as `a:b`).  If c is missing or 0, or more
generally equal to a, then the value is calculated as: `c=1.3 *
max(time)`, where `time` is the hdf5 attribute: the measurement time
of the experiment. A missing value for b be corrected as
`b=(c-a)/100.0`.
