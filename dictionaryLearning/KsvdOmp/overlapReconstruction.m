function [IOut] = overlapReconstruction(patches, patchSize, sz)
patches = double(patches);
row = sz(1);
col = sz(2);
count = 1;
Weight = zeros(row, col);
IMout = zeros(row , col);
idx = [1:size(patches,2)];
[rows,cols] = ind2sub([row-patchSize+1, col - patchSize+1], idx);
for i  = 1:length(cols)
    col = cols(i); row = rows(i);        
    block =reshape(patches(:,count),[patchSize, patchSize]);
    IMout(row:row+patchSize-1,col:col+patchSize-1)=IMout(row:row+patchSize-1,col:col+patchSize-1)+block;
    Weight(row:row+patchSize-1,col:col+patchSize-1)=Weight(row:row+patchSize-1,col:col+patchSize-1)+ones(patchSize);
    count = count+1;
end;
IOut = (IMout)./(Weight);
end