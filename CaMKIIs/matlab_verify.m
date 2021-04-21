model_g = h5info('CaMKIIs.h5');
output_g = h5info('CaMKIIs_out.h5');

nG=length(output_g.Groups);
log_par=h5read('CaMKIIs.h5','/prior/mu');
par=exp(log_par);
npar=length(par);

for i=1:nG
 g_name = output_g.Groups(i).Name;
 fprintf("-----\nExperiment %i (%s)\n",i,g_name);

 cy = h5read('CaMKIIs_out.h5',strcat(g_name,'/state'));
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
 fprintf("difference in the trajectory, aggregated: %g\nand per state variable:\n",rel_err(Y,cy'));
 disp(sum(abs(Y-cy'))./sum(1e-8+abs(Y)));
 ny=length(y0);
 nt=length(t);
 mjac=zeros(ny,ny,nt);
 mjacp=zeros(ny,np,nt);
 cSnorm_t=zeros(1,nt);
 
 for j=1:nt
  mjac(:,:,j)=Jy(t(j),Y(j,:)');
  mjacp(:,:,j)=Jp(t(j),Y(j,:)');
  cSnorm_t(j) = norm(cS(:,:,j));
 end%for
 fprintf("diff between C_jacobian and matlab_jacobian: %g\n",rel_err(mjac,permute(cjac,[2,1,3])));
 fprintf("diff between C_p_jacobian and matlab_p_jacobian: %g\n",rel_err(mjacp,permute(cjacp,[2,1,3])));

 % make a second simulation, with slightly different parameters
 h=1e-6;
 H = randn(np,1)*h;
 p_default=p;
 p=p_default.*H;
 delta_p=p-p_default;
 f=@(t,y) CaMKIIs_vf(t,y,p);
 Jy=@(t,y) CaMKIIs_jac(t,y,p);
 Jp=@(t,y) CaMKIIs_jacp(t,y,p);
 odeset('jacobian',Jy);
 odeset('RelTol',1e-5);
 odeset('AbsTol',1e-6);
 odeset('BDF',true);
 %tspan=[min(t) max(t)];
 [T2,Y2]=ode15s(f,t,y0);
 predicted_Y=Y;
 for j=1:nt
   predicted_Y(j,:) = predicted_Y(j,:) + permute(cS(:,:,j)*delta_p,[3,1,2]);
 end%for
 S_err=rel_err(Y2,predicted_Y);
 fprintf("sensitivity error estimate based on trajectory prediction: %g\n",S_err);
 plot(t,cSnorm_t);
 xlabel('t');
 ylabel('norm(S)');
 %figure(i);
 %plot(T,Y);
 %xlabel('t');
 %ylabel('state');
end%for
