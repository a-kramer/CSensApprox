#!/usr/bin/octave-cli -q
#set(0,"defaultaxeslinewidth",2);
#set(0,"defaultlinelinewidth",2);
set(0,"defaulttextinterpreter","none");
set(0,"defaultaxesfontname","Fira Sans");
set(0,"defaulttextfontname","Fira Sans");

# the actual system is here:
load LinearSystem.txt
assert(logical(exist('A','var')))

load LinearSystem_out.h5
assert(logical(exist('trajectory','var')) && isstruct(trajectory)); # it was loaded from LinearSystem_out.h5


function [x,S]=analytical_solution(A,x0,t,b,u)
  ## Usage: [x,S]=analytical_solution(A,x0,t,b,u)
  ##
  ## Exact solution to the IVP of the affine ode system: A*x + b.^2 + u
  ## with initial value x0, parameters b and input u
  ##
  ## x is a matrix of size length(t) × length(x0)
  ## just like the lsode return value for x
  ##
  ## NOTE: we assume t0 = 0.0
  ##
  nt=length(t);
  nx=length(x0);
  nb=nx;
  x=NA(nt,nx);
  S=NA(nx,nb,nt);
  b2=b.^2;
  c=A\(b2+u);
  B=diag(b);
  C=A\B;
  for i=1:nt
    eA=expm(A*t(i));
    x(i,:)=reshape(eA*(x0+c) - c,1,nx);
    S(:,:,i)=2*(eA*C - C);
  endfor
endfunction

function [x_dp] = parameter_shifted(x,t,S,p,dp)
  dx=S*dp;
  x_dp = x+dx;
endfunction

function [NormMt,varargout]=norm_time_series(M)
  nt=size(M,3);
  NormMt=NA(nt,1);
  for i=1:nt
    NormMt(i)=norm(M(:,:,i));
  endfor
  if (nargout>1)
    rc=NA(nt,1);
    for i=1:nt
      rc(i)=rcond(M(:,:,i));
    endfor
    varargout{1}=rc;
  endif
endfunction

function [fi_t,sv,v,w]=fisher_information(S)
  ## Usage: [fi_t]=fisher_informtaion(S)
  ##
  ## S is a 3-dimensional series of sensitivity matrices
  ## S(:,:,i) is the sensitivity at time point i (t(i))
  n=size(S);
  fi_t=NA(n(2),n(2),n(3));
  w=NA(n(1),n(1),n(3));
  v=NA(n(1),n(2),n(3));
  sv=NA(n(1),n(3));
  for i=1:n(3)
    fi_t(:,:,i)=S(:,:,i)'*S(:,:,i);
    [w(:,:,i),s,v(:,:,i)]=svd(S(:,:,i));
    sv(:,i)=diag(s);
  endfor
endfunction

tt=trajectory.time;
cX=trajectory.state';
tf=100;

figure(1); clf;
subplot(2,2,1); cla;
plot(tt,cX);
title("gsl solution in C");
xlabel("t");
ylabel("x(t)");
xlim([0,1]*tf);
subplot(2,2,2); cla;
u=1.0;
[aX,aS]=analytical_solution(A,x0,tt,b,u);
plot(tt,aX);
title("analytical solution");
xlabel("$t$");
ylabel("$x(t)$");
xlim([0,1]*tf);

subplot(2,2,3); cla;
u=1.0;
tt=trajectory.time;
[aX,aS]=analytical_solution(A,x0,tt,b,u);
semilogy(tt,1e-9+sum(abs(aX-cX),2));
title("Difference between gsl and analytical solution");
xlabel("$t$");
ylabel('$\sum_i |x_i(t) - x_i(t;\texttt{gsl})|$');
xlim([0,1]*tf);


subplot(2,2,4); cla;
load LinearSystem.h5
load time.txt
plot(data.trajectory)
title('lsode solution in \emph{GNU Octave}');
xlabel("$t$");
ylabel("$x(t)$");
xlim([0,1]*tf);

set(gca,"fontname","Fira Sans");
set(gcf,"paperunits","centimeters");
set(gcf,"papersize",[20,16]);
set(gcf,"paperposition",[0,0,20,16]);
print("Trajectory.tex","-dpdflatex");

nb=length(b);
tc=5;
ci=find(tt<tc,1,'last');
cS=permute(trajectory.sensitivity,[2,1,3])(:,1:nb,:);
S_err = norm(sum(abs(cS(:,:,1:ci)-aS(:,:,1:ci)),3))/ci;
printf("average sensitivity error up to t=%g (before it diverges): %g\n",tt(ci),S_err);

dt=diff(tt);

figure(2); clf;
subplot(2,2,1); cla;
tspan=[0.5,1];
in_tspan=tt>=tspan(1) & tt<=tspan(2);
[fi_t,sv1,w,v]=fisher_information(cS(:,:,in_tspan));
boxplot(tt(in_tspan)-dt(in_tspan)/7,sv1,
	"face color",[0.6,0.6,1.0],
	"width",0.25,
	"edge color",[0.2,0.2,0.5],
	"median",[0,0,0]);
[fi_t,sv2,w,v]=fisher_information(aS(:,:,in_tspan));
boxplot(tt(in_tspan)+dt(in_tspan)/7,sv2,
	"face color",[1.0,0.8,0.6],
	"width",0.25,
	"edge color",[0.5,0.2,0.1],
	"median",[0.6,0.1,0.1]);
xlim(tspan);
xlabel('$t$');
ylabel('singular values of $S(t)$ \texttt{svd}');
title('gsl odeiv2 solution (blue) analytical (orange)');

subplot(2,2,2); cla;
tspan=[2,4]
in_tspan=tt>=tspan(1) & tt<=tspan(2);
[cFI_t,sv3,w,v]=fisher_information(cS(:,:,in_tspan));
boxplot(tt(in_tspan)-dt(in_tspan)/7,sv3,
	"face color",[0.6,0.6,1.0],
	"width",0.25,
	"edge color",[0.2,0.2,0.5],
       	"median",[0,0,0]);
[aFI_t,sv4,w,v]=fisher_information(aS(:,:,in_tspan));
boxplot(tt(in_tspan)+dt(in_tspan)/7,sv4,
	"face color",[1.0,0.8,0.6],
	"width",0.25,
	"edge color",[0.5,0.2,0.1],
       	"median",[0.6,0.1,0.1]);
xlim(tspan);
xlabel('$t$');
ylabel('singular values of $S(t)$ \texttt{svd}');
title('gsl odeiv2 solution (blue) analytical (orange)');

[aN,aRC]=norm_time_series(aFI_t);
[cN,cRC]=norm_time_series(cFI_t);
subplot(2,2,3); cla;
plot(tt(in_tspan),aN,";analytical;",tt(in_tspan),cN,";gsl;");
xlabel("$t$");
ylabel('$\|\texttt{FI}(t)\|$');
title("Fisher Information (norm)");

subplot(2,2,4); cla;
l=-round(log10(max(aRC)));
plot(tt(in_tspan),aRC*10^l,";analytical;",tt(in_tspan),cRC*10^l,";gsl;");
xlabel("$t$");
ylabel(sprintf('$\\text{rcond}(\\texttt{FI}) / 10^{-%i}$',l));
title("Fisher Information (reciprocal condition number)");

set(gca,"fontname","Fira Sans");
set(gcf,"paperunits","centimeters");
set(gcf,"papersize",[20,16]);
set(gcf,"paperposition",[0,0,20,16]);
print("FisherInformation.tex","-dpdflatex");
