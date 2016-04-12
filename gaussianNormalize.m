function [ nv ] = gaussianNormalize( v )
%GAUSSIANNORMALIZE Summary of this function goes here
%   Detailed explanation goes here

v(find(isnan(v))) = 0;
nv = bsxfun(@rdivide, (bsxfun(@minus, v , mean(v))), std(v));

end

