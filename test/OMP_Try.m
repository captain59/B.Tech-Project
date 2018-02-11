function S=my_omp(D,X,L,E)
% Orthogonal Matching Pursuit (OMP)
%%%%%%%%   Input:  D: Dictionary (matrix)
%%%%%%%%           X: Signal 
%%%%%%%%           L: The maximum number of non-zero coefficients for each sparse signal
%%%%%%%%           E: Minimum allowed reconstruction error for patch
%%%%%%%%   Output: S: Coefficient vector for sparse representation

[N,K]=size(D); %% N: dim of signal(patch size) ; K: Atoms in dictionary
[N,P]=size(X); %% N: dim of signal(patch size) ; P: Number of patches in each image

if (N~=size(X,1))
    error('Dimension not matched');
end

S=zeros(K,P);
for k=1:1:P
    x=X(:,k); %% k-th patch of image 
    residual=x; %% Residual of X
    Rec_x=zeros(size(x));
    Omega=[];
    
    for j=1:1:L        
        corr=D'*residual; %% Compute correlation between residual and columns of D

        [~, ind]=max(abs(corr)); %% Choose the maximum correlation coefficient index for corresponding columns of D
 
        Omega=union(Omega,ind);
        
        S(Omega,k)=pinv(D(:,Omega))*x;
        
        residual=residual-(D(:,Omega))*S(Omega,k); %% Update the residual
        
%         Rec_x=Rec_x+(D(:,Omega))*S_temp(Omega,k); %% Reconstructed patch
        
    if sum(residual.^2)< E
        break;
    end
    
    end
end