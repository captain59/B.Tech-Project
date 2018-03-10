clc;
close all;
fileName = 'image.png';
Y1 = double(imread(fileName));
mask = double(imread('numMask.png')); % known portionwhite
mask = increaseMaskSize(mask);
Xint = double(zeros(size(Y1)));
lamda = 3.5; learningRate = 0.01; itteration = 500;
for iter=1:itteration
    prior = prior_gmrf(Xint);
    O_est = mask/255;
    %O_est = imcomplement(O_est);
    Y_est = O_est.*Xint;
    Ydiff = Y1 - Y_est;
    gradient = 0.5*-2*O_est.*Ydiff + 5*ssimDerivativeGaussian(Y1, Y_est)+lamda*prior;
    %gradient = 0.5*-2*O_est.*Ydiff+0.5*ssimDerivativeGaussian(Y1, Y_est) + lamda*prior;
    Xint = Xint-learningRate*gradient;
    disp(['Itteratiob ', num2str(iter), ' MSE: ', num2str(mse(Y1, Xint)), ' SSIM: ', num2str(ssim(Y1, Xint))]);
end
figure;
imshow(Xint, []);
figure;
imshow(Y1, []);