function [S] = SparseCodedVector(D, Y, R, delta)
meanY = mean(Y);
w=8;
h=8;
SSIMValues = [];
sparseCodedVectors = [];
rhoValues = meanY;
% for i=1:R
%     rhoValues = [rhoValues, meanY-i*delta, meanY+i*delta];
% end
oneMat = ones(w*h,1);
H = eye(w*h)-(oneMat*oneMat')/(w);
capK = (D'*H*D)/(w*h);%
smallK = (D'*H*Y)/(w*h);
varY = var(Y);
meanDic = (D'*oneMat)/(w*h);
C2 = floor((0.03*255)^2);
for rho=rhoValues                                        
    tou = 0.2;
    Utou = 1.0;
    Ltou = tou;
    epsilon = 0.05;
    count = 0;
    while true
        X = vectorX(capK, smallK, meanDic, tou, rho);
        Stou = tou*(varY+X'*capK*X+C2)-(2*smallK'*X+C2);
        Dtou = Utou-Ltou;
        if Stou>=0 && Dtou<epsilon
            break;
        elseif tou>0.99
            break
        elseif Stou>=0 && Dtou>=epsilon
            tou = (Utou+Ltou)/2;
            Utou = tou;
        else
            tou = (Utou+Ltou)/2;
            Ltou = tou;
        end
        count=count+1;
        if count> 10
            break;
        end
    end
%     disp(count);
    X =  vectorX(capK, smallK, meanDic, tou, rho);
    Ycap = D*X;
    SSIMValues = [SSIMValues, ssim(Ycap, Y)];
    sparseCodedVectors = [sparseCodedVectors, X];
end
[~, ind] = max(SSIMValues);
S = sparseCodedVectors(:, ind);
end