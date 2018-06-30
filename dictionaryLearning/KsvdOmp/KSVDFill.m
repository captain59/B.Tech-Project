clc;
clear all;
close all;
Err_threshold=0.8;
MaxIter=10   ;
Patchsize=[8 8]; %% Patch size
P1=8;
% Input image
[FileName,PathName]=uigetfile('*.jpg','Select image file');
img =imread(fullfile(PathName,FileName));

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
patches = im2col(Vec, Patchsize, 'distinct');
% sample = datasample(1:length(selectPatches), 100000);
% patches = selectPatches(:, sample);
% select patches without occlusions
Blk_Vec = [];
patchLength = length(patches);
patchOccludedIndices = [];
for i=1:patchLength
    wzero = find(~patches(:,i));
    if isempty(wzero)
        Blk_Vec  = [Blk_Vec, patches(:,i)];
    else
        patchOccludedIndices = [patchOccludedIndices, i];
    end
end
disp('Patch extraction complete');
l=1;
ssimValMax = 0;
for Numiter=1:1:MaxIter
    Sparmat = OMP_Try(D,Blk_Vec,r/2, 0.9); %% Sparse matrix optimization    

    learningParameter = 0.1;
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
end
for i=patchOccludedIndices
    targetPatch = patches(:,i);
    wzero = find(~targetPatch);% finds n.o of zero entries
    wz = wzero;
    while(isempty(wzero)==0)
        SparMat = OMP_Try(D, targetPatch, r/2, 0.7);
        targetPatch = D*SparMat;
        wzero = find(~targetPatch);
    end
    patches(wz,i) = abs(targetPatch(wz));
end
Vec_Out = col2im(patches, [8 8], [row col], 'distinct');
figure;
subplot(1,2,1);
imshow(Vec, []);
title('Original Image')
%%% Display reconstructed image
subplot(1,2,2)
imshow(Vec_Out, []);