function [out] = extractMask(img)
[row, col] = size(img);
out = 255*ones(size(img));
for i=1:row
    for j=1:col
        if img(i,j) <= 5
            out(i,j) = 0;
        end
    end
end
end