function [u] = SB_ATV_inpaint(g, mu)
g = double(g);
%g = g(1:end, 1:end-93);
%g = g(1:50, 1:50);
O_est = extractMask(g)/255;
[M, N] = size(g);
[B, Bt, BtB] = DiffOper(sqrt(M*N));
b = zeros(2*M*N,1);
d = b;
u = ones(size(g));
k = 1;
%tol = 1e-3;
lambda = 3;
num = 200;
while k < 10
    fprintf('it. %g ',k);
    %up = u;
    learningParameter = 0.005;
    %error = [];
    for i=1:num
        Ydiff = g - u;
        prior = sparsit(lambda, u, b, d, M, N);
        grad = -2*O_est.*Ydiff + lambda*prior;
        u = u - learningParameter*grad;
        error(i) = sum(sum(Ydiff))/(M*N);
        disp(['Itteration: ', num2str(i)]);
        subplot(1, 3, 1), imshow(g, []); title('Initial Obs');
        subplot(1, 3, 2), imshow(u, []); title('u updated');
        subplot(1, 3, 3), plot(error), title('MSE ');
        drawnow;
    end
    Bub = B*u(:)+b;
    d = max(abs(Bub)-mu/lambda,0).*sign(Bub);
    b = Bub-d;
    %err = norm(up-u)/norm(u);
    %fprintf('err=%g \n',err);
    k = k+1;
end
%fprintf('Stopped because norm(up-u)/norm(u) <= tol=%.1e\n',tol);
end

function [B, Bt, BtB] = DiffOper(N)
D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,1) = [];
D(1,1) = 0;
B = [ kron(speye(N),D) ; kron(D,speye(N)) ];
Bt = B';
BtB = Bt*B;
end

function [priorresult] = sparsit(lambda, u, b, d, M, N)
[B, Bt, BtB] = DiffOper(sqrt(M*N));
u = u(:);
priorresult = lambda*reshape((BtB*u - Bt*(b-d)), M, N);
end

function [u] = gradientDescent(u, g, O_est, lambda, b, d, num)
[M, N] = size(u);
learningParameter = 0.001;
for i=1:num
    Ydiff = g - u;
    prior = sparsit(lambda, u, b, d, M, N);
    grad = -2*O_est.*Ydiff + lambda*prior;
    u = u - learningParameter*grad;
end
end