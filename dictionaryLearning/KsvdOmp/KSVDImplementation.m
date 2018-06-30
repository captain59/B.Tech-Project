clc;
clear all;
close all;


Err_threshold=0.8;
MaxIter=2;
Patchsize=[8 8]; %% Patch size
P1=8;
% Input image
[FileName,PathName]=uigetfile('*.jpg','Select image file');
img =imread(fullfile(PathName,FileName));
% img = imresize(img, 0.5);
if ~isa(img, 'double')
    img=double(img);
end
if size(img,3)>1
    img=rgb2gray(img);
end

Vec = img;
[row,col,band]=size(Vec);

r=ceil(col/2);

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
disp('Initial Dictionary obtained');
% Patches
selectPatches = im2col(Vec, Patchsize);
sample = datasample(1:length(selectPatches), 100000);
patches = selectPatches(:, sample);
% select patches without occlusions
Blk_Vec = [];
patchLength = length(patches);
for i=1:patchLength
    wzero = find(~patches(:,i));
    if isempty(wzero)
        if isempty(Blk_Vec)
            Blk_Vec = patches(:,i);
            disp('Occluded Patch');
        else
            Blk_Vec = cat(2, Blk_Vec, patches(:,i));
        end
    end
end
disp('Patch extraction complete');
l=1;
ssimValMax = 0;
for Numiter=1:1:MaxIter
    Sparmat=OMP_Try(D,Blk_Vec,r/2, 0.9); %% Sparse matrix optimization    

    learningParameter = 0.01;
    for d=1:r
        w_k=find(Sparmat(d,:));
        y = Blk_Vec(:,d);
        atom = D(:, d);
        if ~isempty(w_k)
            for i=w_k
                derivative = Sparmat(d,i)*SSIMDerivative(y, Sparmat(d,i)*D(:,d));
                atom = atom + learningParameter*derivative;
            end
        end
        D(:, d) = atom;
    end   
    disp('finished Trainning dictionary');
    
%     %%% Reconstruction
%     Sparmat=OMPTest(D, patches, r/2, 0.6);
%     Rec_blkmat=D*Sparmat;
%     mse_error=mse(patches ,Rec_blkmat)/band;
%     
%     Vec_Out=col2im(Rec_blkmat,[8 8],[row col], 'distinct');
%     PSNROut=20*log10(255/sqrt(mean((Vec_Out(:)-Vec(:)).^2)));
%     SSIMR=ssim(Vec_Out,Vec);
%     if(SSIMR > ssimValMax)
%         bestImg = Vec_Out;
%     end
%     
%     disp(['Iter: ' num2str(Numiter),'  done,  mse: ', num2str(mse_error),'  psnr: ', num2str(PSNROut),' ssim: ', num2str(SSIMR)]);
%     
%     if SSIMR == 1
%         Numiter=1+MaxIter;
%     end
end
Sparmat=OMPTest(D, selectPatches, r/2, 0.6);
Rec_blkmat=D*Sparmat;
bestImg = overlapReconstruction( Rec_blkmat, 8, size(Vec));
figure;
subplot(1,2,1);
imshow(Vec, []);
title('Original Image')
%%% Display reconstructed image
subplot(1,2,2)
imshow(bestImg, []);
title(strcat(['Retrieved: ',num2str(PSNROut),'dB']));