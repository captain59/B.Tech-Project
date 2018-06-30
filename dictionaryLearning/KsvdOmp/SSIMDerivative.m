function [derivative] = SSIMDerivative(X, Y)
%Taking the SSIM derivative wrt Y
%
Np = length(X);
meanX = mean(X);
meanY = mean(Y);
covX = cov(X);
covY = cov(Y);
covXY = cov(X,Y);
covXY = covXY(1,2);
C1 = floor((0.01*255)^2);
C2 = floor((0.03*255)^2);
A1 = 2*meanX*meanY + C1;
A2 = 2*covXY + C2;
B1 = meanX^2 + meanY^2 + C1;
B2 = covX^2 + covY^2 + C2;
derivative = 2*(A1*B1*(B2*X-A2*Y)+B1*B2*(A2-A1)*meanX+A1*A2*(B1-B2)*meanY)/(Np*B1^2*B2^2);
end