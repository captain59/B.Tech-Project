function [ret] = derivativeSSIM(X, Y)
derivative = [];
Np = length(X);
meanX = mean(X);
meanY = mean(Y);
cov_XY = cov(X,Y);
covXY = cov_XY(1,2);
varX = cov_XY(1,1);
varY = cov_XY(2,2);
C1 = floor((0.01*255)^2);
C2 = floor((0.03*255)^2);
for i=1:Np
    dev = 2*SSIMCalc(X(i), Y(i))*(meanX/(2*meanX*meanY+C1)+(X(i)-meanX)/(2*covXY+C2)-(Y(i)-meanY)/(varX+varY+C2)-2*meanY/(meanX^2+meanY^2+C1));
    derivative = [derivative dev];
end
ret = derivative';
end