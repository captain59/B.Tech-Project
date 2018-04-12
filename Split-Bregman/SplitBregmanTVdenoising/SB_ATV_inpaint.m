g = imread('Lena512.png');
g = double(g);
mu = 20;
%g = g(1:end, 1:end-93);
%g = g(1:50, 1:50);
mask = imread('LenaMask.png');
mask = double(mask);
O_est = mask/255;
g = O_est.*g;
[M, N] = size(g);
[B, Bt, BtB] = DiffOper(sqrt(M*N));
b = zeros(2*M*N,1);
d = b;
u = ones(size(g));
k = 1;
%tol = 1e-3;
lambda = 3;
num = 200;
while k < 30
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
        if error(i) < 0
            learningParameter = 0.95*learningParameter;
        end
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


