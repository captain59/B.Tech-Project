function [S] = ssimDerivativeGaussian(A, ref)
if isempty(A)
    disp('Empty Image');
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
sigma = 2.5;
radius = 2*ceil(3*sigma)+1;
C1 = floor((0.01*255)^2);
C2 = floor((0.03*255)^2);
% Gaussian weighting function
% gaussFilt = getGaussianWeightingFilter(radius,ndims(A));
gaussFilt = fspecial('gaussian', [radius radius] , sigma);

% Weighted-mean and weighted-variance computations
mux = imfilter(A, gaussFilt,'conv','replicate');
muy = imfilter(ref, gaussFilt,'conv','replicate');

muxy = mux.*muy;
mux2 = mux.^2;
muy2 = muy.^2;

sigmax2 = imfilter(A.^2,gaussFilt,'conv','replicate') - mux2;
sigmay2 = imfilter(ref.^2,gaussFilt,'conv','replicate') - muy2;
sigmaxy = imfilter(A.*ref,gaussFilt,'conv','replicate') - muxy;

l = (2*mux.*muy + C1)./(mux2 + muy2 + C1);
cs = (2*sigmaxy + C2)./(sigmax2 + sigmay2 + C2);

%num = (2*muxy + C1).*(2*sigmaxy + C2);
%den = (mux2+ muy2 + C1).*(sigmax2 + sigmay2 +C2);
%val = num./den;
w = fspecial('gaussian', size(A), sigma);
dl = imfilter(2*((muy - mux.*l)./(mux2 + muy2 + C1)), gaussFilt, 'same');
dcs = imfilter(2./(sigmax2 +sigmay2 + C2).*((ref - muy) - cs.*(A - mux)), gaussFilt, 'same');

%meanVal = l.*cs;
%disp(mean(meanVal(:)));
%disp(mean(val(:)));
%disp(ssim(A, ref));

S = -1*(dl.*cs + l.*dcs);
end