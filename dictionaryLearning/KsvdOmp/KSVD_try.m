clc;
clear all;
close all;

Err_threshold=0.000005;
MaxIter=10;
Patchsize=[8 8]; %% Patch size
P1=8;

% Input image
[FileName,PathName]=uigetfile('*.jpg','Select image file');
Vec=imread(fullfile(PathName,FileName));
[row,col,band]=size(Vec);

if ~isa(Vec, 'double')
    Vec=double(Vec);
end
if size(Vec,3)>1
    Vec=rgb2gray(Vec);
end
Blk_Vec=im2col(Vec,Patchsize,'distinct');

% Dictionary size
% r=round(2.0*prod(patch_size));
r=input('Enter the number of atoms in Dictionary: ');
% d_type=input('Type dictionary matrix: \n 1. Uniform distributed entries in range -1 to 1 \n 2. Uniform distributed entries in range 0 to 1 \n 3. Gaussian distributed entries with zeros mean \n');
% if d_type==1
%     D=2*rand(size(Blk_Vec,1),r)-1;
% elseif d_type==2
%     D=rand(size(Blk_Vec,1),r);
% elseif d_type==3
%     D=randn(size(Blk_Vec,1),r);
% end
% D=D.*(ones(size(D,1),1)*(1./sqrt(sum(D.*D)))); %% Normalize dictionary
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
D=cat(2,D,D,D,D,D,D);

l=1;

% Name=sprintf('Iter_%d Sigma_%d Errthresh_%d Denoiseblksize_%d Patchsize_%d',MaxIter,Sigma,Err_threshold,denoise_blksize,8*8);
for Numiter=1:1:MaxIter
    Sparmat=OMP_Try(D,Blk_Vec,r/2,Err_threshold); %% Sparse matrix optimization
    Sparmat_up=zeros(size(Sparmat));
    D_up=zeros(size(D));
    
    if Numiter==1
        Rec_blkmat_init=D*Sparmat;
        mse_error_init=mse(Blk_Vec,Rec_blkmat_init)/band;
        disp(['Initial: ' 0,'   done, initial mse: ', num2str(mse_error_init)]);
    end
    
    for k=1:1:r
        Temp_E=zeros(size(Blk_Vec));
        w_k=find(Sparmat(k,:));
        
        if ~isempty(w_k)
            Sparmat_temp=Sparmat(:,w_k);
            Sparmat_temp(k,:)=0; %% The coeffitients of the element we now improve are not relevant.
            E_k_R=(Blk_Vec(:,w_k) - D*Sparmat_temp); %% Vector of errors that we want to minimize with the new element
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
    
    Vec_Out=col2im(Rec_blkmat,[8 8],[row col],'distinct');
    PSNROut=20*log10(255/sqrt(mean((Vec_Out(:)-Vec(:)).^2)));
    SSIMR=ssim(Vec_Out,Vec);
    
    
    disp(['Iter: ' num2str(Numiter),'  done,  mse: ', num2str(mse_error),'  psnr: ', num2str(PSNROut),' ssim: ', num2str(SSIMR)]);
    
    if SSIMR == 1
        Numiter=1+MaxIter;
    end
end

%%% Display reconstructed image
figure;
subplot(1,2,1)
imshow(Vec,[])
title('Original Image')
subplot(1,2,2)
imshow(Vec_Out,[])
title(strcat(['Noisy: ',num2str(PSNROut),'dB']));