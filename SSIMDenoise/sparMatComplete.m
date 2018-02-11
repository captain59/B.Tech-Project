function [S] = sparMatComplete(D, Y, T, Range, delta)
% T- N.o of sparse coded vector
[~, P] = size(Y);
S = [];
for i=1:P
    Sret = optimalVector2(D, Y(:,i), Range, delta, T);
    S = [S, Sret];
end
end