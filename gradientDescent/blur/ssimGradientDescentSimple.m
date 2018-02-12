function [original, reconstructed, mseMat, ssimMat] = ssimGradientDescentSimple(A, varargin)
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
    learningParameter = 2;
end

X = 100*ones(size(Y));
ssimValMax = -1;
bestImg = X;
disp(['Initial SSIM: ', num2str(ssim(A, X))]);
mseMat = mse(A, X);
ssimMat = ssim(A, X);
for i=1:itterations
    derivative = SSIMDerivative(Y, conv2(X, G, 'same'));
    error = learningParameter*conv2(derivative, G,'same');
    X = X - error;
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