#!/usr/bin/octave-cli -q


A=randn(n=6);
A=round(-A'*A*10)/10;
b=round(rand(n,1)*10)/10;
u=1;

f=@(x,t) A*x+b.^2 + u; # x'=f(x,t)
Jx=@(x,t) A;           # df/dx
Jb=@(x,t) diag(b);     # df/db

x0=ones(n,1);
save -text InitialValue.txt x0

X=lsode({f,Jx},x0,t=linspace(0,100,nt=2^7))

plot(t,X);
xlabel("t");
ylabel("x(t)");
title("state variables");
save time.txt t

save trajectory.txt X
fid=fopen("trajectory.double","w");
fwrite(fid,X,"double");
fclose(fid);

fid=fopen("mu.double","w");
fwrite(fid,log(b),"double");
fclose(fid);


save LinearSystem.txt A b f Jx Jb x0 u

vf=fopen("LinearSystem.vf","w");
fprintf(vf,'<?xml version="1.0" ?>\n');
fprintf(vf,'<VectorField Name="LinearSystem" Description="A*x+b.^2 +u">\n');
for i=1:n
  for j=i:n
    fprintf(vf,'<Constant Name="a%i%i" Value="%g"/>\n',i,j,A(i,j));
  endfor
endfor
for i=1:n
  fprintf(vf,'<Parameter Name="b%i" DefaultValue="%g"/>\n',i,b(i));
endfor
fprintf(vf,'<Parameter Name="u" DefaultValue="%g"/>\n',u);
for i=1:n
  fprintf(vf,'<StateVariable Name="x%i" DefaultInitialCondition="%g" Formula="',i,x0(i));
  for j=1:n
    fprintf(vf,'%sa%i%i*x%i',merge(j==1,"","+"),min(i,j),max(i,j),j);
  endfor
  fprintf(vf,'+b%i*b%i + u',i,i);
  fprintf(vf,'"/>\n');
endfor
fprintf(vf,'</VectorField>\n');
fclose(vf);



