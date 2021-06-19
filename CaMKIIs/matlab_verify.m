% CaMKIIs contains the inputs, initial conditions,
%         used when running the solver
% CaMKIIs_out contains the solutions
% CaMKIIs.h5 -> gsl_odeiv2 solver -> CaMKIIs_out.h5
model_g = h5info('CaMKIIs.h5');
output_g = h5info('CaMKIIs_out.h5');

% Experiments are organized in groups (in the output file)
nG=length(output_g.Groups);
log_par=h5read('CaMKIIs.h5','/prior/mu');
par=exp(log_par);
npar=length(par);

for i=1:6:nG
 g_name = output_g.Groups(i).Name;
 fprintf("-----\nExperiment %i (%s)\n",i,g_name);
 Major=h5readatt('CaMKIIs.h5',strcat('/data',g_name),'major');
 Minor=h5readatt('CaMKIIs.h5',strcat('/data',g_name),'minor');
 cy = h5read('CaMKIIs_out.h5',strcat(g_name,'/state'));
 Status = h5read('CaMKIIs_out.h5',strcat(g_name,'/status'));
 cjac = h5read('CaMKIIs_out.h5',strcat(g_name,'/jac'));
 cjacp = h5read('CaMKIIs_out.h5',strcat(g_name,'/jacp'));
 cS = h5read('CaMKIIs_out.h5',strcat(g_name,'/sensitivity'));
 PHIf = h5read('CaMKIIs_out.h5',strcat(g_name,'/transition_matrix_forward'));
 PHIb = h5read('CaMKIIs_out.h5',strcat(g_name,'/transition_matrix_backward'));
 
 cS=permute(cS,[2,1,3]);
 t=h5read('CaMKIIs_out.h5',strcat(g_name,'/time'));
 u=h5readatt('CaMKIIs.h5',strcat('/data',g_name),'input');
 p=cat(1,par,u);
 np=length(p);
 y0=h5readatt('CaMKIIs.h5',strcat('/data',g_name),'InitialValue');
 f=@(t,y) CaMKIIs_vf(t,y,p);
 Jy=@(t,y) CaMKIIs_jac(t,y,p);
 Jp=@(t,y) CaMKIIs_jacp(t,y,p);
 odeset('jacobian',Jy);
 odeset('RelTol',1e-5);
 odeset('AbsTol',1e-6);
 odeset('BDF',true);
 %tspan=[min(t) max(t)];
 [T,Y]=ode15s(f,t,y0);
 fprintf("difference in the trajectory, aggregated: %g\nand per state variable:\n",norm(mean(rel_err(Y,cy'))));
 % relative error:
 disp(sum(abs(Y-cy'))./sum(1e-8+abs(Y)));
 ny=length(y0);
 nt=length(t);
 mjac=zeros(ny,ny,nt);
 mjacp=zeros(ny,np,nt);
 cSnorm_t=zeros(1,nt);
 
 for j=1:nt
  mjac(:,:,j)=Jy(t(j),Y(j,:)');
  mjacp(:,:,j)=Jp(t(j),Y(j,:)');
 end
 % the norm of the sensitivity is easier to plot than the matrix itself
 for j=1:nt
  cSnorm_t(j) = norm(cS(:,:,j));
 end%for
 fprintf("diff between C_jacobian and matlab_jacobian: %g\n",norm(rel_err(mjac,permute(cjac,[2,1,3]))));
 fprintf("diff between C_p_jacobian and matlab_p_jacobian: %g\n",norm(rel_err(mjacp,permute(cjacp,[2,1,3]))));

 % make a second simulation, with slightly different parameters
 S_err=linearization_error(p,y0,t,cS,Y);
 %S_err(abs(Y)<1e-6)=NaN;
 figure(i); clf;
 semilogy(t,S_err./(1e-6+Y));
 xlabel('time');
 ylabel('Delta_y(t;p)/(eps+y(t;p))');
 title(sprintf('linearization error of E %i.%i',Major,Minor));
 plot_to_file(sprintf('LinErr_%i_%i.png',Major,Minor));
 
 
 l=find(Status<0,1);
 fprintf("average sensitivity error estimate based on trajectory prediction: %g\n",mean(reshape(S_err,1,[]))/mean(reshape(Y,1,[])));
 figure(nG+i); clf;
 subplot(2,1,1);
 plot(t,cSnorm_t);
 if ~isempty(l)
  hold on;
  yL=ylim();
  plot([t(l),t(l)],yL);
 end
 xlabel('t');
 ylabel('norm(S)');
 title(strcat(sprintf('Exp %i.%i with input = [ ',Major,Minor),sprintf('%g ',u),'];'))
 subplot(2,1,2);
 stairs(t,Status);
 xlabel('t');
 ylabel('status');
 plot_to_file(sprintf('SnsvtyNorm_%i_%i.png',Major,Minor));
 %figure(i);
 %plot(T,Y);
 %xlabel('t');
 %ylabel('state');
end%for


function plot_to_file(filename,varargin)
  if (nargin>1)
    sz=varargin{1};
    set(gcf,'papersize',sz);
    set(gcf,'paperposition',cat(2,zeros(1,2),sz));
  end%if
  set(gca,'fontname','Fira Sans');
  set(gca,'fontsize',16);
  saveas(gcf,filename)
end%function
