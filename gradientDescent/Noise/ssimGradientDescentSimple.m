function [reconstructed, original, mseMat, ssimMat] = ssimGradientDescentSimple(A, varargin)
[~, ~, channels] = size(A);
if channels>1
    A = rgb2gray(A);
end
A = double(A);
sigma = 5;

Y = A + sigma*randn(size(A));

if nargin>=2
    itterations = varargin{1};
else
    itterations = 1000;
end
if nargin>=3
    learningParameter = varargin{2};
else
    learningParameter = 2;
end

X = 100*ones(size(Y));
ssimValMax = -1;
bestImg = X;
disp(['Initial SSIM: ', num2str(ssim(A, X))]);
mseMat = mse(A, X);
ssimMat = ssim(A, X);
for i=1:itterations
    derivative = learningParameter*SSIMDerivative(Y, X);
    X = X - derivative;
    ssimVal = ssim(A, X);
    ssimMat = [ssimMat ssimVal];
    mseMat = [mseMat mse(A, X)];
    if ssimVal > ssimValMax
        ssimValMax = ssimVal;
        bestImg = X;
    end
    disp(['Itteration: ', num2str(i),' SSIM: ', num2str(ssimVal)]);
end
reconstructed = bestImg;
original = Y;
imwrite(uint8(reconstructed), 'SSIM_Reconstruction_Simple.jpg');
end