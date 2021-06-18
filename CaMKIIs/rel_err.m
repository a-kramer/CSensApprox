function [E]=rel_err(A,B)
 D=abs(A-B);
 S=max(reshape(A+B,1,[]));
 E=norm(sum(D,3));
 tol=1e-8*(1+max(reshape(S,1,[])));
 E=E/(tol+S);
end%function