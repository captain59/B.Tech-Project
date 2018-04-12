% Master Program to Execute all the files
function masterRunning(img, varargin)
if nargin>=2
    itterations = varargin{1};
else
    itterations = 1000;
end
if nargin>=3
    learningRate = varargin{2};
else
    learningRate = 1;
end
[original, ~, mseMatMSE, ssimMatMSE] = l2gradientdescent(img, itterations, 0.005*learningRate);
[~, ~, mseMatSSIMGaussian, ssimMatSSIMGaussian] = ssimGradientDescent(img, itterations, 5*learningRate);
[~, ~, mseMatSSIMSimple, ssimMatSSIMSimple] = ssimGradientDescentSimple(img, itterations, 8*learningRate);
%Saving Noisy Image
imwrite(uint8(original), 'Original_Noisy.jpg');
% plotting
plotting(mseMatMSE, mseMatSSIMGaussian, mseMatSSIMSimple, 'Variation of MSE with Itterations', 'MSE Error');
hold on
plotting(ssimMatMSE, ssimMatSSIMGaussian, ssimMatSSIMSimple, 'Variation of SSIM with Itterations', 'SSIM');
end