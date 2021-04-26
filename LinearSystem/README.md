# Affine System

This example is using the analytical solution to linear ordinary
differential equations, to solve a (random) affine system. The model
with independent variable _t_ shall have _n_ state variables _x_ and
the same number of parameters, as well as one scalar input _u_.

```bash
##   ẋ = f(x,t;p)
  b[i] = p[i]*p[i]     # i in 0:n-1
     u = p[n]          # 
 dx/dt = A*x + b + u,  # with x(0) = x0
```
where `p` is this model's parameter vector (same size as `x`). As
mentioned in the other examples many solvers have a single slot for
parameters in the expected function interface, so _p_ contains both
the parameters that we presume as being subject to later optimization
_b_ as well as an input to the model _u_, which is a known parameter
that distibguishes different experiments (experimental conditions).

For this example we use [GNU
Octave](https://www.gnu.org/software/octave/index) and command line
tools to create a model and simulation instruction file. And _GNU
Octave_ again to evaluate the simulation results.

## HDF5 Attributes

[GNU Octave](https://www.gnu.org/software/octave/index) has built-in
support for hdf5 files. The function `load()` can read generic hdf5
files and `save()` can write _GNU Octave_ specific files. It's not
very flexible and we need to work around some limitations.

The standard `load()` function does not support reading ATTRIBUTES of
hdf5 data sets. Writing them is not possible either. Instead, `load()`
will import all datasets into the main scope (as nested structure
arrays if necessary, see [hdf5
groups](https://support.hdfgroup.org/HDF5/doc1.6/UG/09_Groups.html),
groups are similar to file system folders/directories).

We solve this problem by using hdf5 command line tools instead when
writing to hdf5 files. Here is an example (in bash):
```bash
H5F=LinearSystem.h5
DATASET=/data/trajectory
h5import trajectory.double -d 6,$((2**7)) -p $DATASET -s 64 -o $H5F
```
However, the hdf5 toolchain does not include a way to write
ATTRIBUTES. h5import only imports DATASETS (as far as we know). 

For this reason we provide the program `h5attr`, it has a similar interface:
```bash
h5attr -d $DATASET -a index -s 0 $H5F
h5attr -d $DATASET -a time time.txt $H5F
```
Here, the value of the attribute can be specified either using the`-s`
option, or by providing an additional file that doesn't end in
`.h5`. A file that ends in `.txt` is assumed to be an ASCII text file,
any other file is assumed to be a binary file which contains double
precision floating point values.

| short |             long | value                 | meaning                  |
|------:|-----------------:|:----------------------|:-------------------------|
|  `-d` |     `--data-set` | string                | hdf5 dataset name        |
|  `-a` |    `--attribute` | string                | name of new attribute    |
|  `-s` | `--string-value` | a string with numbers | attribute value (vector) |
|       |                  | `${FILE}.h5`          | the file to modify       |
|       |                  | `${FILE}.txt`         | attribute value file     |
|       |                  | `${FILE}.${OTHER}`    | attribute value, binary  |

A binary file will be read using `fread()` in C (`stdio.h`). A text
file will be read using `getline()` (GNU extension to gcc). it does
not matter whether the numbers are one value per line, or all in one
line, the attribute value will be a normal one dimensional array
regardless.


### Writing appropriate files in GNU Octave

A text file that is acceptable to `h5attr` can be written using the `save()`
function, by hand or any other method (like `printf()`). A binary file
can be written like this:
```matlab
fid=fopen("trajectory.double","w");
fwrite(fid,X,"double");
fclose(fid);
```

## Verify

The script [verify.m](./verify.m) will compare the analytical solution
of the model to the numerical (gsl/C) calculated trajectories. It will
also make a graphic comparing the analytical and estimated
sensitivity.
