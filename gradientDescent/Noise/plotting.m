function plotting(a, b, c, titleName, YLabel)
N = length(a);
x = 1:N;
if length(a)~= length(b) || length(b)~= length(c)
    disp('Lengths should be equal');
    return
end
fig = figure;
plot(x, a, 'r', x, b, 'g--', x, c, 'b:');
title(titleName);
xlabel('Itterations');
ylabel(YLabel);
legend('MSE', 'SSIM Gaussian', 'SSIM without Gaussian');
saveas(fig, titleName,'png');
end