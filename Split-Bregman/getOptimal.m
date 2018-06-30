function [out] = getOptimal(Iorig, Xint, mask)
maskOne = mask/255;
out = Iorig.*maskOne + Xint.*imcomplement(maskOne);
end