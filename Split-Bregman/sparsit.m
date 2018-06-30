function [priorresult] = sparsit(lambda, u, b, d, M, N)
[B, Bt, BtB] = DiffOper(sqrt(M*N));
u = u(:);
priorresult = lambda*reshape((BtB*u - Bt*(b-d)), M, N);
end