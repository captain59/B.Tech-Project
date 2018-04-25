clc;
close all;
fileName = 'peppers.png';
Y1 = imread(fileName);
Y1 = double(Y1);
mask = imread('peppermask.png');
mask = double(mask);
O_est = mask/255;
Y1 = O_est.*Y1;
%Xinit = double(ones(size(Y1)));
Xinit = Y1;
[M, N] = size(Y1);
%% Initilaization
gammaInit = 5000; gammaTarget = 10; k = 0.9;
lambda = 1; learningRate = 0.05;
Xout = Xinit;
while gammaInit ~= gammaTarget
    % gradient descent
    Xinit = Xout;
    error = [];
    for iter=1:800
        prior = DA_MRF(Xinit, gammaInit);
        Y_est = O_est.*Xinit;
        Ydiff = Y1 - Y_est;
        %gradient = -2*O_est.*(Ydiff)+ lambda*prior;
        % Trying to minimize 1- ssim thus the negative sign is multiplied
        ssimgrad = O_est.*SSIMDerivative(Y1, Y_est);
        gradient = ssimgrad + lambda*prior;
        Xinit = Xinit - learningRate*gradient;
        error(iter) = sum(sum(Ydiff))/(M*N);
        disp(['Itteration ', num2str(iter), ' MSE: ', num2str(mse(Y1, Xinit)), ' SSIM: ', num2str(ssim(Y1, Xinit))]);
        subplot(2, 2, 1), imshow(Y1, []); title('Initial Obs');
        subplot(2, 2, 2), imshow(Xinit, []); title('Xint updated');
        subplot(2, 2, 3), imshow(Y_est, []); title('Y_est updated');
        subplot(2, 2, 4), plot(error); title('Error');
        drawnow;
    end
    gammaInit = max(k*gammaInit, gammaTarget);
    if gammaInit < 100
        learningRate = 0.5*learningRate;
    end
    disp(gammaInit);
    Xout = Xinit;
end
figure;
imshow(Xinit, []);
figure;
imshow(Y1, []);