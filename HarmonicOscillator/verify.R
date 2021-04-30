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
    message(sprintf("## %s",Label))
    message(sprintf("Problem size/sensitivity dimensions: `parameter × state × time`:\n```R"))
    print(dim(cS))
    message("```")

    k <- D[["p"]]
    names(k) <- 'k'

    sv0 <- D[["sv0"]]
    u <- D[["input"]]
    names(u) <- c('c','F')
    names(sv0) <- c("v","y")

    message(sprintf("### Default Values:\n```R"))
    print(k)
    print(sv0)
    message("```")

    w0=sqrt(k)
    sol <- .exact.solution(x,sv0,k,u)
    y.exact <- sol[['y']]
    L <- gsub("[ ]","_",Label)
    
    if (abs(u['F'])<1e-2){
        ## Show that the numerical solution is correct within tolerances
        ## dev.new()
        pdf(sprintf("%s_gsl_vs_exact_solution.pdf",L))
        y.error <- .solution.error(x,sv,sol,k,cS,Label,u)
        dev.off()
        ## Show that the sensitivity can be used to approximate
        ## trajectories with slightly changed parameters
        ## dev.new()
        pdf(sprintf("%s_linearization_error.pdf",L))
        par(mfcol = c(2, 1))
        .linear.sensitivity.error(x,sv,sv0,cS,k,dk=5e-2,u,Label)
        .linear.sensitivity.error(x,sv,sv0,cS,k,dk=1e-1,u,Label)
        dev.off()
        ## We can actually show the sensitivity here:
        pdf(sprintf("%s_sensitivity_gsl_vs_exact_solution.pdf",L))
        y.max <- max(abs(cS[1,2,]))
        plot(x,cS[1,2,],xlab='time',ylab='sensitivity',main='dy(t;p)/dp approximation',ylim=c(-1,1)*y.max)
        lines(x,sol[['S']],lty=1)
        lines(x,sol[['Sfd']],lty=2)
        points(x,sol[['Scif']],lty=3,pch=4)
        legend('topleft',c('approximation','exact','finite differences','Cauchy'),lty=c(NA,1,2,NA),pch=c(1,NA,NA,4))
        dev.off()
    } else {
        ##dev.new()
        pdf(sprintf("%s_compared_to_deSolve.pdf",L))
        .compare.with.deSolve(x,sv,sv0,cS,k,dk=4e-2,u,Label)
        dev.off()
    }
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
         sub="[C] gsl odeiv solution",
         main=Label,
         lty="dashed",
         xlim=c(0,xf),
         ylim=c(-1,1))
    lines(x,sv[2,])
    legend("topright",c("v","y"),lwd=c(1,1),lty=c("dashed","solid"))

    y.exact <- sol[['y']]
    S.exact <- sol[['S']]
    message(sprintf("Analytical sensitivity of the 4th to 8th time point:\n```R"))
    print(S.exact[4:8])
    message("```")
    message(sprintf("Estimated sensitivity of the 4th to 8th time point:\n```R"))
    print(cS[1,2,4:8])
    message("```")

    y.error <- sum(abs(y.exact - sv[2,]))/length(x)
    S.error <- sum(abs(S.exact - cS[1,2,]))/length(x)
    message("The error is calculated as the sum of absolute differences, normalized by the number of
time-points.")
    message(sprintf("The trajectory error from `gsl odeiv2` is %g.",y.error))
    message(sprintf("The approximate sensitivity error is %g.",S.error))
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

.compare.with.deSolve <- function(x,sv,sv0=c(v=0.0,y=1.0),cS,k=exp(-1),dk=1e-2,u=c(0.0,0.0),Label){
        
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
        message(sprintf("### deSolve vs Sensitivity shifted Trajectory:\n Sum of absolute differences per time-point: %g",Diff))
 
        par(mfcol = c(2, 1))
        tm <- sol[, "time"]
        plot(tm, sol[, "v"],
             type = "l",
             xlab = "t",
             ylab = "v",
             main=Label,
             sub=sprintf("deSolve solution v(t;k+dk) with dk=%g",dk))
        lines(x,sv.est[1,],lty=2)
        legend("topright",c("deSolve","v(t;k)+Sv(t;k)*dk"),lty=c(1,2))

        plot(tm, sol[, "y"],
             type = "l",
             xlab = "t",
             ylab = "y",
             main=Label,
             sub=sprintf("average missmatch: %g",Diff))
        lines(x,sv.est[2,],lty=2)
        legend("topright",c("deSolve","y(t;k)+Sy(t;k)*dk"),lty=c(1,2))
    }
}

.linear.sensitivity.error <- function(x,sv,sv0,cS,k,dk,u,Label="linearization error with F=0"){
    dsv <- cS[1,,]*dk
    sv.est <- sv+dsv
    y0 <- sv0[2]
    sol.k0 <- .exact.solution(x,sv0,k,u)
    sol.k <- .exact.solution(x,sv0,k+dk,u)
    xf <- max(x)*1.5
    y.err <- sum(abs(sol.k[['y']]-sv.est[2,]))/length(x)
    plot(x,sol.k[['y']],
         ylab='y',
         type='p',
         pch=1,
         sub=sprintf("average missmatch: %g",y.err),
         main=sprintf("%s at dk=%g",Label,dk),
         xlim=c(0,xf))
    lines(x,sv.est[2,],lty=2)
    w0 <- sqrt(k)
    lines(x,sol.k0[['y']],lty=1)
    lines(x,sol.k0[['y']]+sol.k0[['S']]*dk,lty=3)
    legend("bottomright",
           c("exact y(t;k+dk,c,F=0)",sprintf("y(t;k,c=%g,F=%g)+Sy(t;k,c,F)*dk",u[1],u[2]),"exact y(t;k,c,F=0)"),
           lty=c(0,2,1),
           pch=c(1,NA,NA))
    message("scores: average absolute error of y(t;k+dk) - (y(t;k)+S(t;k)*dk)")
    score=c(finite.differences=sum(abs(sol.k[['y']] - (sol.k0[['y']] + sol.k0[['Sfd']]*dk)))/length(x),
            approximation=sum(abs(sol.k[['y']] - (sv.est[2,])))/length(x),
            exact=sum(abs(sol.k[['y']] - (sol.k0[['y']]+sol.k0[['S']]*dk)))/length(x),
            Cauchy=sum(abs(sol.k[['y']] - (sol.k0[['y']]+sol.k0[['Scif']]*dk)))/length(x))
    print(score)
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
    ##
    ## y(0) = + a*cos(f) = y0
    ## v(0) = - a*cos(f)*r*w - a*sin(f)*sqrt(1-r*r)*w
    ##      = - y0*r*w - y0*tan(f)*sqrt(1-r*r)*w
    ## (v0 + y0*r*w)/(y0*sqrt(1-r*r)*w) = - tan(f)
    ## d/dk ((v0 + y0*r*sqrt(k))/(y0*sqrt(1-r*r)*sqrt(k))) = - dtan(f)/df * df/dk
    ## d/dk (v0/(y0*sqrt(1-r*r)*sqrt(k))) = -2/(cos(2*f+1)) * dfdk
    ## (v0/(y0*sqrt(1-r*r))*(0.5*k^-1.5) = 2/(cos(2*f+1)) * dfdk
    ## => df/dk = (v0*(cos(2*f+1))/(y0*sqrt(1-r*r)*4*sqrt(k)^3) 
    c <- u[1]
    F <- u[2]
    y0 <- sv0['y']
    v0 <- sv0['v']
    ##
    w <- function(k) sqrt(k)
    r <- function(k) c/(2*w(k))        # damping ratio
    srr1 <- function(k) sqrt(1-r(k)^2) # convenience
    ## srr1*w
    srr1w <- function(k) sqrt(k-0.25*c^2)
    ## dw/dk
    dwdk <- function(k) 1/(2*w(k))
    f <- function(k) atan(-(2*v0+y0*c)/(2*y0*srr1w(k)))  # phase
    dfdk <- function(k) cos(2*f(k) + 1)*(2*v0 + y0*c)/(8*y0*srr1w(k)^3)
    a <- function(k) y0/cos(f(k))                       # y0 if v0==0 && r==0
    dadk <- function(k) y0 * (tan(f(k))/cos(f(k))) * dfdk(k)
    y <- function(x,k) a(k)*exp(-0.5*c*x)*cos(srr1w(k)*x + f(k))
    #S <- function(x,k) exp(-0.5*c*x) * (dadk(k)*cos(srr1w(k)*x + f(k)) - a(k)*sin(srr1w(k)*x + f(k))*(x/(2*srr1w(k)) + dfdk(k)))
    S <- function(x,k) y(x,k)*(tan(f(k))*dfdk(k) - tan(srr1w(k)*x+f(k))*(0.5*x/srr1w(k) + dfdk(k)))
    h <- 1e-6
    Sfd <- (-y(x,k+2*h)+8*y(x,k+h)-8*y(x,k-h)+y(x,k-2*h))/(12*h)
    ## Cauchy's integral formula
    R <- 3
    g <- function(s,k) (k+R*(cos(s)+sin(s)*1i))
    n <- length(x)
    I <- vector(mode='numeric',length=n)
    N <- 361
    s.0.2pi <- seq(0,2*pi,length.out=N)
    ds <- mean(diff(s.0.2pi))
    for (i in 1:n){
        xi <- x[i]
        G <- function(s,k) y(xi,g(s,k))/(g(s,k)-k)^2 * ((g(s,k)-k)*1i)
        I[i] <- sum(G(s.0.2pi,k)[1:N-1])*ds
    }
    Cauchy  <- Re(I/(2*pi*1i))
    return(list(y=y(x,k),S=S(x,k),Sfd=Sfd,Scif=Cauchy))
}

load_and_verify <- function(ModelName="HarmonicOscillator"){
    output.h5 <- sprintf("%s_out.h5",ModelName)
    h5f <- h5file(output.h5)
    for (g in names(h5f)){
        gl <- gsub("([A-Z])"," \\L\\1\\E",g,perl=TRUE)
        Default <- .defaults(sprintf("%s.h5",ModelName),g)
        ## print(Default)
        h5g <- openGroup(h5f,g)
        verify(h5g,Default,gl)
        h5close(h5g)
    }
    h5close(h5f)
}
