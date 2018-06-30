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
        index = i;
    end
end
Dchange = D(:, index);
%initial sparse coded vector
Stemp(index) = 1;
% finding the T vectors
for i=1:(T-1)
    arr = 1:K;
    L = ~ismember(arr, index);
    arr = arr.*L;
    iter = find(arr);
    ind = 1;
    ssimMaxVal = -1;
    for j=iter
        Dfunc = [Dchange, D(:,j)];
        Sret = SparseCodedVector(Dfunc, Y, R, delta, [Stemp, 1]');
        stemp = SSIMCalc(Dfunc*Sret, Y);
        if stemp > ssimMaxVal
            ssimMaxVal = stemp;
            ind = j;
            SretMax = Sret;
        end
    end
    Dchange = [Dchange, D(:, ind)];
    index = [index, ind];
    Stemp = SretMax';
end
Sfin = SparseCodedVector(Dchange, Y, R, delta, Stemp');
S = zeros(K,1);
S(index) = Sfin;
end