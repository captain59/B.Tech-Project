% DA_MRF prior
function [ prior ] = DA_MRF(X, gamma)
[rows, cols] = size(X);
smooth_arr = zeros(size(X));
for i=2:rows-1
    for j=2:cols-1        
        smooth_arr(i,j)=2*(X(i,j)-X(i,j-1))*exp((-(X(i,j)-X(i,j-1))^2)/gamma) + ...
        2*(X(i,j)-X(i,j+1))*exp((-(X(i,j)-X(i,j+1))^2)/gamma) + ...
        2*(X(i,j)-X(i-1,j))*exp((-(X(i,j)-X(i-1,j))^2)/gamma)+ ... 
        2*(X(i,j)-X(i+1,j))*exp((-(X(i,j)-X(i+1,j))^2)/gamma);
    end
end
prior = smooth_arr;
end