%gmrf/smooth prior
function [ gmrf ] = prior_gmrf( X )
[rows, cols] = size(X);
smooth_arr = zeros(size(X));
for i=2:rows-1
    for j=2:cols-1        
        smooth_arr(i,j)=2*(X(i,j)-X(i,j-1))+2*(X(i,j)-X(i,j+1))...
            +2*(X(i,j)-X(i-1,j))+2*(X(i,j)-X(i+1,j));
    end
end
gmrf = smooth_arr;
end