function [S] = OMP_Try(D, Y, L, Upper_Threshold, Lower_Threshold)
%       D - the dictionary (its columns MUST be normalized).
%       Y - the signals to represent
%       L - the max. number of coefficients for each signal.
% output arguments: 
%       A - sparse coefficient matrix.
[N P] = size(Y);
[N K] = size(D);

if (N~=size(Y, 1))
    error('Dimensions not matched');
end
S = zeros(K, P);
for k=1:P
    y = Y(:,k);
    residual = y;
    %Rec_y = zeros(size(y));
    Omega = [];
    j=0;
    ssimVal = (Upper_Threshold+Lower_Threshold)/2;
    while(ssimVal < Upper_Threshold && ssimVal > Lower_Threshold  && j<L)
        corr = D'*residual;
        [~, ind] = max(abs(corr));
        Omega = union(Omega, ind);
        S(Omega, k) = pinv(D(:, Omega))*y;
        residual = y - D(:, Omega)*S(Omega, k);
        ssimVal = ssim(D(:, Omega)*S(Omega, k), y);
        j=j+1;
    end 
    %disp([num2str(ssimVal),' ',num2str(k)]);
end
return