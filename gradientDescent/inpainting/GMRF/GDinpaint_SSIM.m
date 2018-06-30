%% Does Inpainting using SSIM as the loss functions
%% Uses GMRF as the method for inpainting
clc;
close all;
fileName = 'cloudSample.png';
Y1 = imread(fileName);
Y1 = double(Y1);
load('mask.mat');
mask = mask2;
%mask = imread('peppermask.png');
%mask = imcomplement(mask);
mask = double(mask);
%Xint = double(ones(size(Y1)));
O_est = mask/255;
Y1 = O_est.*Y1;
[M, N] = size(Y1);
Xinit = Y1;
lambda = 1; learningRate = 0.0005;
Xout = Xinit;
itteration = 2500;
for iter=1:itteration
    prior = prior_gmrf(Xinit); 
    %O_est = imcomplement(O_est);
    Y_est = O_est.*Xinit;
    Ydiff = Y1 - Y_est;
    % never use gaussian
    %gradient = 0.5*-2*O_est.*Ydiff + 10*ssimDerivativeGaussian(Y1, Y_est)+lamda*prior;
    ssimgrad = O_est.*SSIMDerivative(Y1, Y_est);
    gradient = 4500*ssimgrad + lambda*prior;
    Xinit = Xinit - learningRate*gradient;
    error(iter) = sum(sum(Ydiff))/(M*N);
    disp(['Itteration ', num2str(iter), ' MSE: ', num2str(mse(Y1, Xinit)), ' SSIM: ', num2str(ssim(Y1, Xinit))]);
    subplot(2, 2, 1), imshow(Y1, []); title('Initial Obs');
    subplot(2, 2, 2), imshow(Xinit, []); title('Xint updated');
    subplot(2, 2, 3), imshow(Y_est, []); title('Y_est updated');
    subplot(2, 2, 4), plot(error); title('Error');
    drawnow;
    Xout = Xinit;
end
figure;
imshow(Xout, []);
figure;
imshow(Y1, []);