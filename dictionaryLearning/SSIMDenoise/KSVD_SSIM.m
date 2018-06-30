function [D_up] = KSVD_SSIM(Sparmat, Y, D)
[~, r] = size(D);
D_up = zeros(size(D));
for k=1:r
    w_k = find(Sparmat(k, :));
    if ~isempty(w_k)
        SparmatTemp = Sparmat(:, w_k);
        SparmatTemp(k, :) = 0;
        E_k_R = (Y(:,w_k)-D*SparmatTemp);
        [U, S, V] = svds(E_k_R);
        D_up(:,k) = U(:, 1);
    end
end
end