clc;
close all;
fileName = 'cloudSample.png';
Y1 = imread(fileName);
Y1 = double(Y1);
load('mask.mat');
mask = mask2;
%mask = imread('newMask.png');
mask = double(mask);
O_est = mask/255;
Y1 = O_est.*Y1;
%Xinit = double(ones(size(Y1)));
Xinit = Y1;
[M, N] = size(Y1);
%% Initilaization
gammaInit = 10000; gammaTarget = 1; k = 0.8;
lambda = 1; learningRate = 0.05;
Xout = Xinit;
masterpath = 'cloudresults/';
path = ['cloudresults/lr_',num2str(learningRate), '_lambda_', num2str(lambda)];
mkdir(path);
while gammaInit ~= gammaTarget
    % gradient descent
    Xinit = Xout;
    error = [];
    for iter=1:400
        prior = DA_MRF(Xinit, gammaInit);
        Y_est = O_est.*Xinit;
        Ydiff = Y1 - Y_est;
        %gradient = -2*O_est.*(Ydiff)+ lambda*prior;
        mseGrad = -2*O_est.*(Ydiff);
        % Trying to minimize 1- ssim thus the negative sign is multiplied
        ssimgrad = O_est.*SSIMDerivative(Y1, Y_est);
        gradient = 4500*ssimgrad + lambda*prior;
        Xinit = Xinit - learningRate*gradient;
        error(iter) = sum(sum(Ydiff))/(M*N);
        disp(['Itteration ', num2str(iter), ' MSE: ', num2str(mse(Y1, Xinit)), ' SSIM: ', num2str(ssim(Y1, Xinit)), ' Gamma: ', num2str(gammaInit)]);
        subplot(2, 2, 1), imshow(Y1, []); title('Initial Obs');
        subplot(2, 2, 2), imshow(Xinit, []); title('Xint updated');
        subplot(2, 2, 3), imshow(Y_est, []); title('Y_est updated');
        subplot(2, 2, 4), plot(error); title('Error');
        drawnow;
    end
    gammaInit = max(k*gammaInit, gammaTarget);
%     if gammaInit < 100
%         learningRate = 0.9*learningRate;
%     end    
    Xout = Xinit;
    imwrite(uint8(Xout), [path, '/',num2str(gammaInit), '.png']);
end
figure;
imshow(Xout, []);
figure;
imshow(Y1, []);