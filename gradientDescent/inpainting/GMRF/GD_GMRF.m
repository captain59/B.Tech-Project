%%Depth map text inpainting+gradient descent+gmrf prior
% from ground truth disparity value
%12/4/17 Sukla
clc; clear all;  close all; 
tic
%%
Y1=double(imread('image.png'));
% Y1_org = load('masked_adirondack');
 O_est=imread('numMask.png');
  [rows,cols]=size(Y1);
  Xint=zeros(size(Y1));
lamda= 20;  grdnt_step=0.001;  grdnt_iter=1000;

    for iter = 1:grdnt_iter
      %str = strcat('D:\Depth_Inpainting\GMRF\result\iter',...
      %num2str(iter),'.png');
     % imwrite(uint8(255*(Xint/max(Xint(:)))),str)
                   
                    prior=prior_gmrf(Xint);
                    
                    O_est=Y1/255;           %O     
                    %O_est = imcomplement(O_est);
                    Y_est=O_est.*Xint;        %Ox
                    Ydiff=Y1-Y_est;           %y-Ox
                    gradient= -2*O_est.*(Ydiff)+lamda*prior;
                    
                    Xint=Xint - grdnt_step*gradient;   
                                    
                    Ydiff1(iter)=sum(sum(Ydiff.^2))/(rows*cols);
                    

                   % subplot(221),imshow(uint8(X),[]);title('Ground truth ');
                    subplot(222),imshow(uint8(Y1),[]);title('Initial obs ');
                    subplot(223),imshow(uint8(Xint),[]);title('Gmrf depth text inpaint');
                    subplot(224),plot(Ydiff1),title('Error plot ');
                    drawnow;
                    iter
    end
             