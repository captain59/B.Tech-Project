function [orig, S, mseMat, ssimMat] = ssimGradientDescent(A, varargin)
[~, ~, channels] = size(A);
if channels>1
    A = rgb2gray(A);
end
A = double(A);
sigma = 1.5;
s = 2*ceil(3*sigma)+1;
G = fspecial('gaussian', [s s], sigma);
Y = conv2(A, G,'same');

if nargin>=2
    itterations = varargin{1};
else
    itterations = 1000;
end
if nargin>=3
    learningParameter = varargin{2};
else
    learningParameter = 0.5;
end

X = 100*ones(size(Y));
ssimValMax = -1;
bestImg = X;
disp(['Initial SSIM: ', num2str(ssim(Y, X))]);
mseMat = mse(Y, X);
ssimMat = ssim(Y, X);
for i=1:itterations
    derivative = ssimDerivativeGaussian(Y, conv2(X, G,'same'));
    error = learningParameter*conv2(derivative, G,'same');
    X = X + error;
    ssimVal = ssim(Y, X);
    ssimMat = [ssimMat ssimVal];
    mseMat = [mseMat mse(Y, X)];
    if ssimVal > ssimValMax
        ssimValMax = ssimVal;
        bestImg = X;
    end
    disp(['Itteration: ', num2str(i),' SSIM: ', num2str(ssimVal)]);
end
S = bestImg;
orig = Y;
imwrite(uint8(S), 'SSIM_Reconstruction_Gaussian.jpg');
end