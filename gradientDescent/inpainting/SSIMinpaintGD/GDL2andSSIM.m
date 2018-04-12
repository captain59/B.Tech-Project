%%
clc;
close all;
clear all;
%%
Y1 = imread('peppers.png');
Y1 = double(Y1);
mask = imread('peppermask.png');
%mask = imcomplement(mask);
mask = double(mask);
%%
O_est = mask/255;
Y1 = O_est.*Y1;
Xint = double(ones(size(Y1)));
%Xint = Y1;
lamda = 1; learningRate = 0.005; mu = 20; itteration = 2500;
[M, N] = size(Y1);
%%
for iter=1:itteration
    prior = prior_gmrf(Xint); 
    %O_est = imcomplement(O_est);
    Y_est = O_est.*Xint;
    Ydiff = Y1 - Y_est;
    % Never use gaussian
    %gradient = 0.5*-2*O_est.*Ydiff + 10*ssimDerivativeGaussian(Y1, Y_est)+lamda*prior;
    ssimGrad = 0.5*O_est.*SSIMDerivative(Y1, Y_est);
    gradient = -0.5*O_est.*Ydiff + lamda*prior;
    Xint = Xint - learningRate*gradient;
    disp(['Itteration ', num2str(iter), ' MSE: ', num2str(mse(Y1, Xint)), ' SSIM: ', num2str(ssim(Y1, Xint))]);
    subplot(1, 3, 1), imshow(Y1, []); title('Initial Obs');
    subplot(1, 3, 2), imshow(Xint, []); title('Xint updated');
    subplot(1, 3, 3), imshow(Y_est, []); title('Y_est updated');
    drawnow;
end
%%
figure;
imshow(Xint, []);
figure;
imshow(Y_est, []);