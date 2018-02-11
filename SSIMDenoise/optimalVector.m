function [S] = optimalVector(D, Y, R, delta, T)
% T - the n.o of sparse coded vectors
ssimMaxVal = -1;
Dchange = [];
index = [];
[~, K]= size(D);
% initial dictionary choosen
for i=1:K
    stemp = SSIMCalc(D(:,i), Y);
    if stemp>ssimMaxVal
        ssimMaxVal = stemp;
        Dchange = D(:,i);
        index = i;
    end
end
% finding the T vectors
ssimMaxVal = -1;
for i=1:(T-1)
    arr = 1:K;
    L = ~ismember(arr, index);
    arr = arr.*L;
    iter = find(arr);
    ind = 1;
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