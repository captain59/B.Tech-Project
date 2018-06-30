clc;
clear all;
close all;

Err_threshold=0.00005;
MaxIter=100;
Patchsize=[8 8]; %% Patch size
P1=8;

% Input image
[FileName,PathName]=uigetfile(['*.jpg':'*.png'],'Select image file');
img=imread(fullfile(PathName,FileName));
Vec = img;
if ~isa(Vec, 'double')
    Vec=double(Vec);
end
if size(Vec,3)>1
    Vec=rgb2gray(Vec);
end
img = Vec;
Vec = img + 12*randn(size(img));
[row col band] = size(img);
% Patches
patches = im2col(Vec, Patchsize);
% sample = datasample(1:length(selectPatches), 100000);
% patches = selectPatches(:, sample);
% select patches without occlusions
Blk_Vec = patches;
% patchLength = length(patches);
% patchOccludedIndices = [];
% for i=1:patchLength
%     wzero = find(~patches(:,i));
%     if isempty(wzero)
%         Blk_Vec  = [Blk_Vec, patches(:,i)];
%     else
%         patchOccludedIndices = [patchOccludedIndices, i];
%     end
% end
disp('Patch extraction complete');
r = 256;
% Dictionary size
%% Normalize dictionary
Pn=ceil(sqrt(r));
DCT=zeros(P1,Pn);
for k=0:1:Pn-1
    V=cos([0:1:P1-1]'*k*pi/Pn);
    if k>0
        V=V-mean(V);
    end
    DCT(:,k+1)=V/norm(V);
end
D=kron(DCT,DCT);

l=1;
disp(ssim(img, Vec));
for Numiter=1:MaxIter
    Sparmat= OMPTest(D, Blk_Vec,10, 0.8, 0.6);
    Sparmat_up=zeros(size(Sparmat));
    D_up=zeros(size(D));
    disp('Running');
    if Numiter==1
        Rec_blkmat_init=D*Sparmat;
        mse_error_init=mse(Blk_Vec,Rec_blkmat_init)/band;
        disp(['Initial: ' 0,'   done, initial mse: ', num2str(mse_error_init)]);
    end
    
    for k=1:r
        Temp_E=zeros(size(Blk_Vec));
        w_k=find(Sparmat(k,:));
        
        if ~isempty(w_k)
            Sparmat_temp=Sparmat(:,w_k);
            Sparmat_temp(k,:)=0; % The coeffitients of the element we now improve are not relevant.
            E_k_R=(Blk_Vec(:,w_k) - D*Sparmat_temp); % Vector of errors that we want to minimize with the new element
            [U,S,V]=svds(E_k_R);
            D_up(:,k)=U(:,1);
            Sparmat_temp=V(:,1)*S(1,1);
            Sparmat_up(k,w_k)=Sparmat_temp';
        else
            Temp_E=Blk_Vec-D*Sparmat;
            E_k=sum(Temp_E.^2);
            [~,i] = max(E_k);
            D_temp=Blk_Vec(:,i);
            D_temp=D_temp./sqrt(D_temp'*D_temp);
            D_up(:,k)=D_temp.*sign(D_temp(1));
            Sparmat_up(k,:)=0;
        end
    end
    
    D=D_up;
    Sparmat=Sparmat_up;
    
    disp('finished Trainning dictionary');
    
    %%% Reconstruction
    Rec_blkmat=D*Sparmat;
    mse_error=mse(Blk_Vec,Rec_blkmat)/band;
    
    Vec_Out=overlapReconstruction(Rec_blkmat, P1, [row col]);
    PSNROut=20*log10(255/sqrt(mean((Vec_Out(:)-img(:)).^2)));
    SSIMR=ssim(Vec_Out,img);
    
    
    disp(['Iter: ' num2str(Numiter),'  done,  mse: ', num2str(mse_error),'  psnr: ', num2str(PSNROut),' ssim: ', num2str(SSIMR)]);
    
    if SSIMR == 1
        Numiter=1+MaxIter;
    end
end
%% Removing Noise
% for i=patchOccludedIndices
%     targetPatch = patches(:,i);
%     wzero = find(~targetPatch);% finds n.o of zero entries
%     wz = wzero;
%     while(isempty(wzero)==0)
%         SparMat = OMPTest(D,Blk_Vec, 10, 0.7, 0.65);
%         targetPatch = D*SparMat;
%         wzero = find(~targetPatch);
%     end
%     patches(wz,i) = abs(targetPatch(wz));
% end
% Vec_Out = col2im(patches, [8 8], [row col], 'distinct');
% S = OMPerr(D, Blk_Vec, Err_threshold);
% Out = D*S;
% Vec_Out = overlapReconstruction(Out, P1, [row col]);
%% Display reconstructed image
figure;
subplot(1,2,1)
imshow(Vec,[])
title('Noisy Image')
subplot(1,2,2)
imshow(Vec_Out,[])
title('Reconstructed Image');