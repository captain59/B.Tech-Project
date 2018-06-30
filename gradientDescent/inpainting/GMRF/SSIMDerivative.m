function [derivative] = SSIMDerivative(Y, ref)
%Taking the SSIM derivative wrt Y
ref = double(ref);
Y = double(Y);
%Np = length(ref(:));
meanX = mean(ref(:));
meanY = mean(Y(:));
cov_XY = cov(ref,Y);
covXY = cov_XY(1,2);
varX = cov_XY(1,1);
varY = cov_XY(2,2);
C1 = floor((0.01*255)^2);
C2 = floor((0.03*255)^2);
A1 = 2*meanX*meanY + C1;
A2 = 2*covXY + C2;
B1 = meanX^2 + meanY^2 + C1;
B2 = varX + varY + C2;
derivative = -2*(A1*B1*(B2*ref-A2*Y)+B1*B2*(A2-A1)*meanX+A1*A2*(B1-B2)*meanY)/(B1^2*B2^2);
end