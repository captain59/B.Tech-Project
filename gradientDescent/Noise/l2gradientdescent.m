function [reconstructed, original, mseMat, ssimMat] = l2gradientdescent(A, varargin)
[~, ~, channels] = size(A);
if channels>1
    A = rgb2gray(A);
end
A = double(A);
sigma = 5;
s = 2*ceil(3*sigma)+1;

Y = A + sigma*randn(size(A));
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
mseValMax = mse(A, X);
mseMat = mseValMax;
ssimMat = ssim(A, X);
bestImg = X;
for i=1:itterations
    derivative = L2Derivative(Y, X);
    X = X - learningParameter*derivative;
    mseVal = mse(A, X);
    ssimMat = [ssimMat ssim(A, X)];
    mseMat = [mseMat mseVal];
    if mseVal < mseValMax
        mseValMax = mseVal;
        bestImg = X;
    end
    disp(['Itteration: ', num2str(i),' MSE: ', num2str(mseVal)]);
end
reconstructed = bestImg;
original = Y;
imwrite(uint8(reconstructed), 'MSE_reconstruction.jpg');
end

function [S] = L2Derivative(Y, X)
S = -2.0*(Y - X);
end