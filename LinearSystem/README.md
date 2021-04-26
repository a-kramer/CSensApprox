# Affine System

This example is using the analytical solution to linear ordinary differential equations, to solve the affine system:
```
b[i]=p[i]^2
dx/dt = A*x + b + u,  with x(0) = x0
```
where `p` is the parameter vector (same size as `x`).

For this example we use [GNU
Octave](https://www.gnu.org/software/octave/index).

## HDF5 Attributes

GNU Octave has built-in support for hdf5 files. The function `load()`
can read generic hdf5 files and `save()` can write _GNU Octave_
specific files. 

However, `load()` does not support reading ATTRIBUTES of hdf5 data
sets. Writing them is also not possible. Instead, `load()` will import
all datasets (into nested structs if necessary).

We solve this problem by using hdf5 command line tools instead when
writing to hdf5 files. Here is an example (in bash):
```bash
H5F=LinearSystem.h5
DATASET=/data/trajectory
h5import trajectory.double -d 6,$((2**7)) -p $DATASET -s 64 -o $H5F
```

However, the hdf5 toolchain does not include a way to write
ATTRIBUTES. h5import only imports DATASETS. 

For this reason we provide the program `h5attr`, it has a similar interface:
```
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
```
fid=fopen("trajectory.double","w");
fwrite(fid,X,"double");
fclose(fid);
```

## Verify

The script [verify.m](./verify.m) will compare the analytical solution
of the model to the numerical (gsl/C) calculated trajectories. It will
also make a graphic comparing the analytical and estimated
sensitivity.
