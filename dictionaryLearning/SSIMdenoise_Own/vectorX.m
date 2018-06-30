function [S] = vectorX(capK, smallK, meanDic, tou, rho)
denominator = (meanDic'/(capK))*meanDic;
lambda = 2*((smallK'/(capK))*meanDic-tou*rho)/denominator;
S = ((capK)\(2*smallK-lambda*meanDic))/(2*tou);
end