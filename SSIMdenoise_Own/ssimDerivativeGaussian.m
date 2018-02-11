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
sigma = 5;

C1 = floor((0.01*255)^2);
C2 = floor((0.03*255)^2);
% Gaussian weighting function
% gaussFilt = getGaussianWeightingFilter(radius,ndims(A));
gaussFilt = fspecial('gaussian', size(A) , sigma);

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
S = -(dl*cs + l*dcs);
end