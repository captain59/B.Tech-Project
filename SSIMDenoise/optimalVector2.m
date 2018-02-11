function [S] = optimalVector2(D, Y, R, delta, T)
% T - the n.o of sparse coded vectors
Dchange = [];
index = [];
[~, K]= size(D);
% finding the T vectors
ssimMaxVal = -1;
for i=1:T
    iter = 1:K;
    for j=iter
        Dfunc = [Dchange, D(:,j)];
        Sret = SparseCodedVector(Dfunc, Y, R, delta);
        stemp = SSIMCalc(Dfunc*Sret, Y);
        if stemp > ssimMaxVal
            ssimMaxVal = stemp;
            ind = j;
        end
    end
    Dchange = [Dchange, D(:, ind)];
    index = [index, ind];
end
Sfin = SparseCodedVector(Dchange, Y, R, delta);
S = zeros(K,1);
S(index) = Sfin;
end