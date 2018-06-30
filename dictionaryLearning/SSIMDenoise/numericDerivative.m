function [ret] = derivativeSSIM(X, Y)
epsilon = 0.0005;
Np = length(X);
derivative = [];
for i=1:Np
    dev = (ssim(X,Y+epsilon)-ssim(X, Y-epsilon))/(2*epsilon);
    derivative = [derivative dev];
end
ret = derivative';
end