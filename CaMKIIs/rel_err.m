function [E]=rel_err(A,B)
 D=abs(A-B);
 S=abs(A+B);
 E=norm(sum(D,3));
 tol=1e-8*(1+max(reshape(S,1,[])));
 E=E/(tol+norm(sum(S,3)));
end%function