% Function to Calculate the M-SSIM derivative
function [S] = MSSIMDerivative(A, ref, sigmaVector)
% sigmaVector is a vector containing the different sigma values used
% M is the size of the sigmaVector
if isempty(A)
    disp('Empty Image');
    return;
end
if isempty(sigmaVector)
    disp('Sigma Vector is empty');
    return;
end

if isa(A,'int16') % int16 is the only allowed signed-integer type for A and ref.
    % Add offset for signed-integer types to bring values in the
    % non-negative range.
    A = double(A) - double(intmin('int16'));
    ref = double(ref) - double(intmin('int16'));
elseif isinteger(A)
    A = double(A);
    ref = double(ref);
end

M = length(sigmaVector);

C1 = floor((0.01*255)^2);
C2 = floor((0.03*255)^2);
[height, width] = size(A);
dcsVector = zeros(height, width, M);
csVector = zeros(1, M);

for i=1:M
    % Gaussian Weigting Function
    gaussFilt = fspecial('gaussian', size(A) , sigmaVector(i));
    % Weighted-mean and weighted-variance computations
    mux = sum(sum(imfilter(A, gaussFilt,'conv','replicate')));
    muy = sum(sum(imfilter(ref, gaussFilt,'conv','replicate')));

    sigmax2 = sum(sum(imfilter(A.^2,gaussFilt,'conv','replicate'))) - mux.^2;
    sigmay2 = sum(sum(imfilter(ref.^2,gaussFilt,'conv','replicate'))) - muy.^2;
    sigmaxy = sum(sum(imfilter(A.*ref,gaussFilt,'conv','replicate'))) - mux.*muy;

    l = (2*mux*muy + C1)/(mux.^2 + muy.^2 + C1);
    cs = (2*sigmaxy + C2)/(sigmax2 + sigmay2 + C2);
    
    dl = 2*gaussFilt*(muy - mux*l)/(mux.^2 + muy.^2 + C1);
    dcs = 2/(sigmax2 +sigmay2 + C2)*gaussFilt.*((ref - muy) - cs*(A - mux));
    %Append to vector
    dcsVector(:,:,i) = dcs;
    csVector(i) = cs;
end

dMSSIM = dl;
for i=1:M
    dMSSIM = dMSSIM + l*dcsVector(:,:,i)/csVector(i);
end
dMSSIM = dMSSIM*prod(csVector);
S = -1*dMSSIM;
end