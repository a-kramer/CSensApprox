function [E]=rel_err(A,B)
 D=abs(A-B);
 S=mean(abs(A+B)/2,3);
 E=mean(D,3);
 tol=1e-8*(1+max(reshape(S,1,[])));
 E=E./(tol+S);
end%function