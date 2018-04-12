function [u] = ATV_Denoising(y, mu, lambda, num)
dx = 0;
dy = 0;
bx = 0;
by = 0;
u = ones(size(y));
for i=1:num
    [Gx, Gy] = imgradientxy(u, 'intermediate');
    up = gradientDescent(y, u, mu, dx, Gx, bx, dy, Gy, by, lambda, num);
    [Gx, Gy] = imgradientxy(up, 'intermediate');
    dx = shrink(Gx + bx, 1/lambda);
    dy = shrink(Gy + by, 1/lambda);
    bx = bx + Gx - dx;
    by = by + Gy - dy;
    disp([ 'Itteration: ',num2str(i), ' MSE: ', num2str(mse(y, up))]);
    u = up;
end
end

function [u] = gradientDescent(y, u, mu, dx, Gx, bx, dy, Gy, by, lambda, num)
for i=1:num
    grad = -mu*(y - u) + -lambda*(dx - Gx - bx) + -lambda*(dy -Gy - by);
    % update u
    u = u - 0.001*lambda*grad;
end
end

function [ret] = shrink(x, lam)
    ret = sign(x).*max(abs(x) - lam, 0);
end