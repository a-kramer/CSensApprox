#!/usr/bin/Rscript
library(hdf5r)
h5f <- h5file("HarmonicOscillator.h5",mode="w")
h5g <- createGroup(h5f, "data")
time <- seq(0,16,length.out=32)
z <- time*0

Data <- matrix(c(z,z),2,32)

h5d <- createDataSet(h5g, "NoDampingNoDrivingForce", Data)
h5attr(h5d,"index") <- 0
h5attr(h5d,"time") <- time
h5attr(h5d,"input") <- c(0.0,0.0)
h5attr(h5d,"InitialValue") <- c(0.0,1.0)

h5d <- createDataSet(h5g, "WithDampingButNoDrivingForce", Data)
h5attr(h5d,"index") <- 1
h5attr(h5d,"time") <- time
h5attr(h5d,"input") <- c(1.0,0.0)
h5attr(h5d,"InitialValue") <- c(0.0,1.0)

h5d <- createDataSet(h5g, "WithDampingAndDrivingForce", Data)
h5attr(h5d,"index") <- 2
h5attr(h5d,"time") <- time
h5attr(h5d,"input") <- c(1.0,0.01)
h5attr(h5d,"InitialValue") <- c(0.0,1.0)

LogParameter <- -1 
h5g <- createGroup(h5f, "prior")
h5d <- createDataSet(h5g, "mu", LogParameter)


h5flush(h5f)
h5close(h5f)
