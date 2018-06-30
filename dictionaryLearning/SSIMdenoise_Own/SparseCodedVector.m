function [S] = SparseCodedVector(D, Y, R, delta, Stemp)
w=8;
h=8;
SSIMValues = [];
sparseCodedVectors = [];
oneMat = ones(w*h,1);
meanDic = (D'*oneMat)/(w*h);
try
    meanY = meanDic'*Stemp;
catch
    disp(meanDic');
    disp('Stemp');
    disp(Stemp);
end
rhoValues = meanY;
for i=1:R
    rhoValues = [rhoValues, meanY-i*delta, meanY+i*delta];
end
H = eye(w*h)-(oneMat*oneMat')/(w*h);
if (rcond(H) < 0.001)
    H = H + 0.01*eye(size(H));
end
capK = (D'*(H')*H*D)/(w*h);%
smallK = (D'*(H')*H*Y)/(w*h);
varY = Stemp'*capK*Stemp;
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
    SSIMValues = [SSIMValues, SSIMCalc(Ycap, Y)];
    sparseCodedVectors = [sparseCodedVectors, X];
end
[~, ind] = max(SSIMValues);
S = sparseCodedVectors(:, ind);
end