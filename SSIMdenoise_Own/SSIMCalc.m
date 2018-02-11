function [ssimVal] = SSIMCalc(x, y)
[row , col] = size(x);
x = double(x);
y = double(y);
X = reshape(x, [1 row*col]);
Y = reshape(y, [1 row*col]);
meanX = mean(X);
meanY = mean(Y);
cov_XY = cov(X,Y);
covXY = cov_XY(1,2);
varX = cov_XY(1,1);
varY = cov_XY(2,2);
C1 = floor((0.01*255)^2);
C2 = floor((0.03*255)^2);
A = (2*meanX*meanY+C1);
B = (2*covXY+C2);
C = (meanX^2+meanY^2+C1);
D = (varX+varY+C2);
ssimVal = (A*B)/(C*D);
end