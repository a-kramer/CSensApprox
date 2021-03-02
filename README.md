# CSensApprox
Experiments with Sensitivity Approximation (SA) in C using the gsl ode solvers.
This program will solve an ODE model for parameters stored in an hdf5 file and any input vector stored in the same hdf5 file.

The depenencies should be minimal, for now this is:

1. GNU Scientific Library
2. HDF5 and specifically HDF5LT

The model will be loaded from a shared library, to keep the model subject to change.
