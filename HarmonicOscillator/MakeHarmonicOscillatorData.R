#!/usr/bin/Rscript
library(hdf5r)
h5f <- h5file("HarmonicOscillator.h5",mode="w")
h5g <- createGroup(h5f, "data")
time <- seq(0,16,length.out=32)

LogParameter <- -1 
w=sqrt(exp(LogParameter))
z <- cos(w*time) + rnorm(32,0,0.05)

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

h5d <- createDataSet(h5g, "WithLittleDampingButNoDrivingForce", Data)
h5attr(h5d,"index") <- 2
h5attr(h5d,"time") <- time
h5attr(h5d,"input") <- c(0.3,0.0)
h5attr(h5d,"InitialValue") <- c(0.0,1.0)


h5d <- createDataSet(h5g, "WithDampingAndSmallDrivingForce", Data)
h5attr(h5d,"index") <- 3
h5attr(h5d,"time") <- time
h5attr(h5d,"input") <- c(1.0,0.01)
h5attr(h5d,"InitialValue") <- c(0.0,1.0)

h5d <- createDataSet(h5g, "WithDampingAndDrivingForce", Data)
h5attr(h5d,"index") <- 4
h5attr(h5d,"time") <- time
h5attr(h5d,"input") <- c(0.8,0.3)
h5attr(h5d,"InitialValue") <- c(0.0,1.0)

h5g <- createGroup(h5f, "prior")
h5d <- createDataSet(h5g, "mu", LogParameter)
h5close(h5g)
h5close(h5f)
