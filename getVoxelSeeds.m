function [ seeds ] = getVoxelSeeds( p, d )
%GETVOXELSEEDS Summary of this function goes here
%   Detailed explanation goes here
% p : Nx3 points
% d : double radius
% seeds :Mx3 seeds

minP = min(p);
maxP = max(p);

radius = d;

[Xq,Yq,Zq] = meshgrid(minP(1):radius:maxP(1),minP(2):radius:maxP(2),minP(3):radius:maxP(3));

seeds = [Xq(:) Yq(:) Zq(:)];

end

