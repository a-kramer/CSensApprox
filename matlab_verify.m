model_g = h5info('CaMKIIs.h5')
output_g = h5info('CaMKIIs_out.h5')
addpath("~/Documents/CaMKIIs/matlab");

nG=length(output_g.Groups);
log_par=h5read('CaMKIIs.h5','/prior/mu');
par=exp(log_par);

for i=1:nG
 g_name = output_g.Groups(i).Name;
 cy = h5read('CaMKIIs_out.h5',strcat(g_name,'/state'));
 cjac = h5read('CaMKIIs_out.h5',strcat(g_name,'/jac'));
 cjacp = h5read('CaMKIIs_out.h5',strcat(g_name,'/jacp'));
 t=h5read('CaMKIIs_out.h5',strcat(g_name,'/time'));
 u=h5readatt('CaMKIIs.h5',strcat('/data',g_name),'input');
 p=cat(1,par,u);
 np=length(p);
 y0=h5readatt('CaMKIIs.h5',strcat('/data',g_name),'InitialValue');
 f=@(t,y) CaMKIIs_vf(t,y,p);
 Jy=@(t,y) CaMKIIs_jac(t,y,p);
 Jp=@(t,y) CaMKIIs_jacp(t,y,p);
 odeset('jacobian',Jy);
 %tspan=[min(t) max(t)];
 [T,Y]=ode15s(f,t,y0);
 fprintf("difference in the trajectory: %g\n",norm(Y-cy'));
 %disp(Y-cy');
 ny=length(y0);
 nt=length(t);
 mjac=zeros(ny,ny,nt);
 mjacp=zeros(ny,np,nt);
 for j=1:nt
  mjac(:,:,j)=Jy(t(j),Y(j,:)');
  mjacp(:,:,j)=Jp(t(j),Y(j,:)');
 end%for
 fprintf("diff between C_jacobian and matlab_jacobian: %g\n",norm(sum(mjac-permute(cjac,[2,1,3]),3)));
 fprintf("diff between C_p_jacobian and matlab_p_jacobian: %g\n",norm(sum(mjacp-permute(cjacp,[2,1,3]),3)));
 %figure(i);
 %plot(T,Y);
 %xlabel('t');
 %ylabel('state');
end%for