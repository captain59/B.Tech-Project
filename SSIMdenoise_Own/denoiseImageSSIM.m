clc;
clear;
close all;
Patchsize=[8 8]; %% Patch size
P1=8;
% Input image
[FileName,PathName]=uigetfile('*.jpg','Select image file');
img =imread(fullfile(PathName,FileName));

% img = randi([0 255], [16 16]);
% img = mat2gray(img);
if ~isa(img, 'double')
    img=double(img);
end
if size(img,3)>1
    img=rgb2gray(img);
end
sigma = 30;
Vec = img + sigma*randn(size(img));
disp(['Initial ssim ', num2str(SSIMCalc(Vec, img))]);
[row,col,band]=size(Vec);

r=64;

%initial DCT dictionary
Pn=ceil(sqrt(r));
DCT=zeros(P1,Pn);
for k=0:Pn-1
    V=cos([0:P1-1]'*k*pi/Pn);
    if k>0
        V=V-mean(V);
    end
    DCT(:,k+1)=V/norm(V);
end
D=kron(DCT,DCT);
DCT = D;
disp('Initial Dictionary obtained');

% Patches
Blk_Vec = im2col(Vec, Patchsize);

disp('Patch extraction complete');
Range = 6;
delta = 5;
%Dictionary Update
[~, r] = size(D);
maxIter = 10;
learningParameter = 0.9;
ssimBest = -1;
bestimg = [];
for j=1:maxIter
    disp(['Itteration ', num2str(j)]);
    if j==1
        SparMat = sparMatComplete(D, Blk_Vec, 10, 5, 6);
    else
        SparMat = S;
    end
    for d=1:r
        w_k=find(SparMat(d,:));
        atom = D(:, d);
        if ~isempty(w_k)
            for i=w_k
                derivative = SparMat(d,i)*SSIMDerivative(Blk_Vec(:,i), SparMat(d,i)*D(:,d));
                atom = atom + learningParameter*derivative;
            end
        end
        D(:, d) = atom;
    end   
    disp('finished Trainning dictionary');
    S = sparMatComplete(D, Blk_Vec, 10, Range, delta);
    p = D*S;
    Vec_Out = overlapReconstruction(p, P1, [row col]);
    PSNROut=20*log10(255/sqrt(mean((Vec_Out(:)-Vec(:)).^2)));
    SSIMR=SSIMCalc(Vec_Out,Vec);
    disp(['  psnr: ', num2str(PSNROut),' ssim: ', num2str(SSIMR)]);
    if SSIMR > ssimBest
        bestimg = Vec_Out;
        ssimBest = SSIMR;
    end
end
