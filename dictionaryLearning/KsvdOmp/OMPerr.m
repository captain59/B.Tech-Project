function [A]=OMPerr(D,X,errorGoal)
%=============================================
% Sparse coding of a group of signals based on a given 
% dictionary and specified number of atoms to use. 
% input arguments: D - the dictionary
%                  X - the signals to represent
%                  errorGoal - the maximal allowed representation error for
%                  each siganl.
% output arguments: A - sparse coefficient matrix.
%=============================================
[n,P]=size(X);
[n,K]=size(D);
E2 = errorGoal;
maxNumCoef = n/2;
A = sparse(size(D,2),size(X,2));
for k=1:P,
    x=X(:,k);
    residual=x;
	indx = [];
	a = [];
	j = 0;
    ssimVal = 0;
    while ssimVal < E2 & j < maxNumCoef,
		j = j+1;
        proj=D'*residual;
        pos=find(abs(proj)==max(abs(proj)));
        pos=pos(1);
        indx(j)=pos;
        a=pinv(D(:,indx(1:j)))*x;
        residual=x-D(:,indx(1:j))*a;
		ssimVal = ssim(D(:,indx(1:j))*a, x);
   end;
   if (length(indx)>0)
       A(indx,k)=a;
   end
end;
return;
