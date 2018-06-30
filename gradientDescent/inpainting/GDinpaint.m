clc;
close all;
fileName = 'peppers.png';
Y1 = imread(fileName);
Y1 = double(Y1);
%load('mask.mat');
mask = imread('peppermask.png');
%mask = imcomplement(mask);
mask = double(mask);
Xint = double(ones(size(Y1)));
lamda = 2; learningRate = 5; itteration = 2500;
O_est = mask/255;
Y1 = O_est.*Y1;
[M, N] = size(Y1);
for iter=1:itteration
    prior = prior_gmrf(Xint); 
    %O_est = imcomplement(O_est);
    Y_est = O_est.*Xint;
    Ydiff = Y1 - Y_est;
    % never use gaussian
    %gradient = 0.5*-2*O_est.*Ydiff + 10*ssimDerivativeGaussian(Y1, Y_est)+lamda*prior;
    gradient = learningRate*SSIMDerivative(Y1, Y_est) + 0.005*lamda*prior;
    Xint = Xint - gradient;
    error = sum(sum(Ydiff))/(M*N);
    if error < 0
        learningRate = 0.95*learningRate;
    end
    disp(['Itteration ', num2str(iter), ' MSE: ', num2str(mse(Y1, Xint)), ' SSIM: ', num2str(ssim(Y1, Xint))]);
    subplot(1, 3, 1), imshow(Y1, []); title('Initial Obs');
    subplot(1, 3, 2), imshow(Xint, []); title('Xint updated');
    subplot(1, 3, 3), imshow(Y_est, []); title('Y_est updated');
    drawnow;
end
figure;
imshow(Xint, []);
figure;
imshow(Y1, []);