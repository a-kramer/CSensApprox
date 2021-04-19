#!/usr/bin/Rscript
library(hdf5r)
h5f <- h5file("HarmonicOscillator_out.h5")
h5g <- openGroup(h5f,"NoDampingNoDrivingForce")
y <- readDataSet(h5g[["state"]])
x <- readDataSet(h5g[["time"]])
cS <- readDataSet(h5g[["sensitivity"]])
jac <- readDataSet(h5g[["jac"]])
jacp <- readDataSet(h5g[["jacp"]])

plot(x,y[1,],type="l",xlab="t",ylab="y, v",main="HarmonicOscillator",sub="no damping, no driving force",lty="dashed",ylim=c(-1,1))
lines(x,y[2,])
legend("topright",c("v","y"),lwd=c(1,1),lty=c("dashed","solid"))

dk <- 1e-2;
dy <- cS[1,,]*dk

y.est <- y+dy

source("HarmonicOscillator_demo.R")

D <- t(sol[,"y"]) - y.est[2,]
