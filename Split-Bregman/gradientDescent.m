function [u] = gradientDescent(u, g, O_est, lambda, b, d, num)
[M, N] = size(u);
learningParameter = 0.001;
for i=1:num
    Ydiff = g - u;
    prior = sparsit(lambda, u, b, d, M, N);
    grad = -2*O_est.*Ydiff + lambda*prior;
    u = u - learningParameter*grad;
end
end