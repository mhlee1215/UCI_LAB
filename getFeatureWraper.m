function [ desc1, key1 ] = getFeatureWraper( v1, c1 )
%GETFEATUREWRAPER Summary of this function goes here
%   Detailed explanation goes here

params.normalRadius=0.1;
params.searchRadius=0.1;
params.searchK = 0;

[v1_s,~,c1_s] = uniformSubSample(v1, 100, c1);
v1_s = v1_s';
c1_s = c1_s';

[desc1, key1] = getFeatures(v1_s, c1_s, params);

end

