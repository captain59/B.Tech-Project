function [reconstructed, original, mseMat, ssimMat] = ssimL2gradientdescent(A, varargin)
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
    learningParameter = 0.01;
end
X = 100*ones(size(Y));

disp(['Initial MSE: ', num2str(mse(A, X))]);
mseValMax = mse(Y, X);
mseMat = mseValMax;
ssimMat = ssim(Y, X);
bestImg = X;
for i=1:itterations
    derivative = 0.5*L2Derivative(Y, X, G) + 0.5*conv2(SSIMDerivative(Y, X), G, 'same');
    X = X - learningParameter*derivative;
    mseVal = mse(A, X);
    ssimVal = ssim(A, X);
    ssimMat = [ssimMat ssim(A, X)];
    mseMat = [mseMat mseVal];
    if mseVal < mseValMax
        mseValMax = mseVal;
        bestImg = X;
    end
    disp(['Itteration: ', num2str(i),' MSE: ', num2str(mseVal), ' SSIM: ', num2str(ssimVal)]);
end
reconstructed = bestImg;
original = Y;
imwrite(uint8(reconstructed), 'ssimL2Rec.jpg');
end

function [S] = L2Derivative(Y, X, G)
S = -2.0*conv2((Y - conv2(X, G,'same')), G,'same');
end