#!/usr/bin/Rscript

## NOTE: c() is the concatenation function in R
## NOTE: t() is the transpose function in R
## dots in names have no special meaning in R: in x.a, x is not a structure

library(hdf5r)
source("HarmonicOscillator.R")

verify <- function(h5g,D,Label){
    sv <- readDataSet(h5g[["state"]])
    x <- readDataSet(h5g[["time"]])
    cS <- readDataSet(h5g[["sensitivity"]])
    jac <- readDataSet(h5g[["jac"]])
    jacp <- readDataSet(h5g[["jacp"]])
    message("sensitivity approximation, dimensions (parameter × state × time): ")
    print(dim(cS))
    message("defaults:")
    k <- D[["p"]]
    names(k) <- 'k'
    print(k)
    sv0 <- D[["sv0"]]
    u <- D[["input"]]
    names(u) <- c('c','F')
    names(sv0) <- c("v","y")
    print(sv0)
    w0=sqrt(k)
    sol <- .exact.solution(x,sv0,k,u)
    y.exact <- sol[['y']]
    
    ##dev.new()
    L <- gsub("[ ]","_",Label)
    pdf(sprintf("%s_gsl_vs_exact_solution.pdf",L))
    y.error <- .solution.error(x,sv,sol,k,cS,Label,u)
    dev.off()
    
    ##dev.new()
    pdf(sprintf("%s_linearization_error.pdf",L))
    par(mfcol = c(2, 1))
    .linear.sensitivity.error(x,sv,sv0,cS,k,dk=2e-2,u,Label)
    .linear.sensitivity.error(x,sv,sv0,cS,k,dk=4e-2,u,Label)
    dev.off()
    
    ##dev.new()
    pdf(sprintf("%s_compared_to_deSolve.pdf",L))
    .compare.with.deSolve(x,sv,sv0,cS,k,dk=4e-2,u)
    dev.off()
}

.defaults <- function(DataFile,Name){
    h5f <- h5file(DataFile)
    h5g <- openGroup(h5f,"prior")
    mu <- readDataSet(h5g[["mu"]])
    h5close(h5g)
    h5g <- openGroup(h5f,"data")
    sv0 <- h5attr(h5g[[Name]],"InitialValue")
    u <-  h5attr(h5g[[Name]],"input")
    h5close(h5g)
    h5close(h5f)
    D <- list(p=exp(mu),sv0=sv0,input=u)
    return(D)
}

.solution.error <- function(x,sv,sol,k,cS,Label,u){
    ## analytical solution:
    ## y''=-k*y
    ## y=y0*cos(w*t)
    ## y'=-y0*sin(w*t)*w
    ## y''=-y0*cos(w*t)*w*w = -y*w*w
    ## k=w*w
    ##
    ## dy/dk = -y0*sin(w*t)*t/(2*sqrt(k))
    par(mfcol = c(2, 1))
    F=u[2]
    xf=max(x)*1.2
    plot(x,sv[1,],
         type="l",
         xlab="t",
         ylab="y, v",
         main="[C] gsl odeiv solution",
         sub=Label,
         lty="dashed",
         xlim=c(0,xf),
         ylim=c(-1,1))
    lines(x,sv[2,])
    legend("topright",c("v","y"),lwd=c(1,1),lty=c("dashed","solid"))

    y.exact <- sol[['y']]
    S.exact <- sol[['S']]
    message("exact sensitivity of the 4th to 8th time point:")
    print(S.exact[4:8])

    message("estimated sensitivity of the 4th to 8th time point:")
    print(cS[1,2,4:8])

    y.error <- sum(abs(y.exact - sv[2,]))/length(x)
    S.error <- sum(abs(S.exact - cS[1,2,]))/length(x)
    message(sprintf("[gsl odeiv] the trajectory error is %g",y.error))
    message(sprintf("[sensitivity approximation] the sensitivity error is %g",S.error))
    plot(x,y.exact,
         lty="solid",
         ylim=c(-1,1),
         lwd=1,
         sub=sprintf("analytical solution (missmatch: %g)",y.error),
         main=Label)
    lines(x,sv[2,])
    legend("topright",c("exact for F=0",sprintf("gsl odeiv with F=%g",F)),lty=c(0,1),pch=c(1,NA))
    return(y.error)
}

.compare.with.deSolve <- function(x,sv,sv0=c(v=0.0,y=1.0),cS,k=exp(-1),dk=1e-2,u=c(0.0,0.0)){
        
    dsv <- cS[1,,]*dk
    sv.est <- sv+dsv
    if(require("deSolve")){

        parameters = c(
            k = k+dk,
            c = u[1],
            F = u[2]
        )
        #print(parameters)

        sol = ode(y = sv0, times = x, func = HarmonicOscillator, parms = parameters,
                  jactype = "fullusr", jacfunc = HarmonicOscillator_jac,
                  atol = 1e-4, rtol = 1e-4)
        y.R <- t(sol[,"y"])
        Diff <- sum(abs(y.R - sv.est[2,]))/length(x)
        message(sprintf("[missmatch] deSolve vs sensitivity estimate: %g",Diff))
 
        par(mfcol = c(2, 1))
        tm <- sol[, "time"]
        plot(tm, sol[, "v"], type = "l",
             xlab = "t", ylab = "v",main=sprintf("deSolve solution at dk=%g",dk))
        plot(tm, sol[, "y"], type = "l",
             xlab = "t", ylab = "y",main=sprintf("comparison to deSolve"),sub=sprintf("sum(abs(sol[y](t;k+dk) - y(k)+[dy/dk]*dk)) = %g",Diff))
        lines(x,sv.est[2,],lty=2)
        legend("topright",c("deSolve","y(t;k)+[dy/dk]*dk"),lty=c(1,2))
    }
}

.linear.sensitivity.error <- function(x,sv,sv0,cS,k,dk,u,Label="linearization error"){
    dsv <- cS[1,,]*dk
    sv.est <- sv+dsv
    y0 <- sv0[2]
    sol.k0 <- .exact.solution(x,sv0,k,u)
    sol.k <- .exact.solution(x,sv0,k+dk,u)
    xf <- max(x)*1.2
    y.err <- sum(abs(sol.k[['y']]-sv.est[2,]))/length(x)
    plot(x,sol.k[['y']],
         type='p',
         pch=1,
         sub=sprintf("sum(abs(y(t;k+dk) - (y(t;k) + [dy/dk]*dk)))/nt = %g",y.err),
         main=sprintf("%s at dk=%g",Label,dk),
         xlim=c(0,xf))
    lines(x,sv.est[2,],lty=2)
    w0 <- sqrt(k)
    lines(x,sol.k0[['y']],lty=1)
    legend("bottomright",
           c("y(t;k+dk)","y(t;k)+S_y_k*dk","y(t;k)"),
           lty=c(0,2,1),
           pch=c(1,NA,NA))
}

.exact.solution <- function(x=c(0,1,2,3),sv0,k=1.0,u=c(0.0,0.0)){
    ## analytical solution: x is the independent variable (time)
    ## y''=-k*y - c*y'
    ## 2*r*w = c
    ##
    ## solution:
    ##             y=a*exp(-r*w*x)*cos(sqrt(1-r*r)*w*x+f)
    ##
    ## check:
    ##
    ## y'=-a*exp(-r*w*x)*(r*w)*cos(sqrt(1-r*r)*w*x+f)
    ##    -a*exp(-r*w*x)      *sin(sqrt(1-r*r)*w*x+f) * sqrt(1-r*r)*w
    ## y' = - y*(r*w)
    ##      - a*exp(-r*w*x)*sin(sqrt(1-r*r)*w*x+f)*sqrt(1-r*r)*w
    ##
    ## y'' = - y'*r*w
    ##       + a*exp(-r*w*x)*(r*w)*sin(sqrt(1-r*r)*w*x+f)*sqrt(1-r*r)*w
    ##       - y * (1-r*r)*w*w
    ##
    ## y'' = - y'*r*w
    ##       + (-y*r*r*w*w - y'*r*w)
    ##       - y * (1-r*r)*w*w
    ##
    ## y'' = - 2*y'*r*w
    ##       - y * r*r*w*w
    ##       - y*ww
    ##       + y * r*r*w*w
    ##
    ## y'' = -y'*c - y*k      [success]
    ##                 ^
    ##               (k=w*w)
    ##
    ## dy/dk = -a*sin(w*t+f)*t/(2*sqrt(k))
    ##
    ## y(0) = + a*cos(f) = y0
    ## v(0) = - a*cos(f)*r*w - a*sin(f)*sqrt(1-r*r)*w
    ##      = - y0*r*w - y0*tan(f)*sqrt(1-r*r)*w
    ## (v0 + y0*r*w)/(y0*sqrt(1-r*r)*w) = - tan(f)
    c <- u[1]
    F <- u[2]
    ##
    w <- sqrt(k)
    r <- c/(2*w)                   # damping ratio
    srr1 <- sqrt(1-r*r)            # convenience
    ## dw/dk
    dwdk <- 1/(2*w)
    y0 <- sv0['y']
    v0 <- sv0['v']
    f <- atan(-(v0+y0*r*w)/(y0*srr1*w))  # phase
    a <- y0/cos(f)                      # y0 if v0==0 && r==0

    y <- a*exp(-r*w*x)*cos(srr1*w*x + f)
    S <- -a*exp(-r*w*x)*(r*dwdk*x)*cos(srr1*w*x + f) - a*exp(-r*w*x)*sin(srr1*w*x + f)*(srr1*dwdk*x)
    return(list(y=y,S=S))
}

load_and_verify <- function(ModelName="HarmonicOscillator",GroupName=c("NoDampingNoDrivingForce","WithDampingButNoDrivingForce","WithDampingAndDrivingForce")){
    output.h5 <- sprintf("%s_out.h5",ModelName)
    print(output.h5)
    print(GroupName[1])
    h5f <- h5file(output.h5)

    for (g in GroupName){
        gl <- gsub("([A-Z])"," \\L\\1\\E",g,perl=TRUE)
        Default <- .defaults(sprintf("%s.h5",ModelName),g)
        print(Default)
        h5g <- openGroup(h5f,g)
        verify(h5g,Default,gl)
        h5close(h5g)
    }
    h5close(h5f)
}
