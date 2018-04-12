function [mask] = increaseMaskSize(img)
[row, col] = size(img);
mask = img;
for i=3:row-3
    for j = 3:col-3
        if img(i,j) == 0
            mask(i-2, j) = 0;
            mask(i-1, j) = 0;
            mask(i+1, j) = 0;
            mask(i+2, j) = 0;
        end
    end
end
end