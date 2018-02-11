clc;
clear;
close all;
Patchsize=[8 8]; %% Patch size
P1=8;
arr = randi([1 64], [64 2]);
Y = arr + max(std(arr))*randn(size(arr));
D = randi([1 64],[64 64]);
arr = double(arr);
D = double(D);
[~, r] = size(D);
[~, sz] = size(Y);
maxIter = 30;
learningParameter = 0.9;
for j=1:maxIter
    disp(['Itteration ', num2str(j)]);
    if j==1
        SparMat = sparMatComplete(D, Y, 10, 5, 6);
    else
        SparMat = S;
    end
    for v=1:sz
        for d=1:r
            w_k=find(SparMat(d,:));
            y = Y(:,v);
            atom = D(:, d);
            if ~isempty(w_k)
                for i=w_k
                    derivative = SparMat(d,i)*SSIMDerivative(y, SparMat(d,i)*D(:,d));
                    atom = atom + learningParameter*derivative;
                end
            end
            D(:, d) = atom;
        end   
    end
%     D = KSVD_SSIM(SparMat, Y, D);
    disp('finished Trainning dictionary');
    S = sparMatComplete(D, Y, 10, 5, 6);
    p = D*S;
    MSE = mse(p, arr);
    PSNR = 10*log10(255^2/MSE);
    SSIM = SSIMCalc(p, arr);
    disp(['Itteration ', num2str(j),' psnr: ', num2str(PSNR),' ssim: ', num2str(SSIM)]);
end