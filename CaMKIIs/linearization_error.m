function [S_err]=linearization_error(p,y0,t,cS,y)
  h=1e-6;
  np=length(p);
  nt=length(t);
  ny=length(y0);
  H = randn(np,1)*h;
  p_default=p;
  p=p_default+H;
  p(55:59)=p_default(55:59); % input has to stay the same
  delta_p=p-p_default;
  f=@(t,y) CaMKIIs_vf(t,y,p);
  Jy=@(t,y) CaMKIIs_jac(t,y,p);
  Jp=@(t,y) CaMKIIs_jacp(t,y,p);
  odeset('jacobian',Jy);
  odeset('RelTol',1e-5);
  odeset('AbsTol',1e-6);
  odeset('BDF',true);
  [T,Y]=ode15s(f,t,y0);
  predicted_y=NaN(nt,ny);
  for j=1:nt
    predicted_y(j,:) = y(j,:) + permute(cS(:,:,j)*delta_p,[3,1,2]);
  end%for
  S_err=abs(Y-predicted_y);
end%function
