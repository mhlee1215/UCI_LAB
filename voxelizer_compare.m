function [ dist, seeds, overlap, d ] = voxelizer_compare( p1, c1, p2, c2, params )
%VOXELIZER_COMPARE Summary of this function goes here
%   Detailed explanation goes here


if ~isfield(params, 'density')
    density = 30;
else
    density = params.density;
end

if ~isfield(params, 'compAll')
    compAll = 0;
else
    compAll = params.compAll;
end

p_all = [p1 ; p2];
d = max(max(p_all)-min(p_all)) / density;

if ~isfield(params, 'seeds')
    seeds = getVoxelSeeds( p_all, d);
else
    seeds = params.seeds;
end

v_params.seeds = seeds;
v_params.d = d;
if ~isfield(params, 'useFilter')
    v_params.useFilter = 0;
else
    v_params.useFilter = params.useFilter;
end
v1 = voxelizer( p1, c1, v_params );
v2 = voxelizer( p2, c2, v_params );

comp_valid = (max(v1')' > 0).*(max(v2')' > 0);
overlap = sum(comp_valid);
diff = v1-v2;
if compAll == 1
    %do nothing
else
    diff = diff(find(comp_valid));
end
dist = sum(sum(diff'*diff));

end

