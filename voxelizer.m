function [ voxelized, seeds ] = voxelizer( p, c, n, params )
% %VOXELIZER Summary of this function goes here
% %   Detailed explanation goes here


seeds = params.seeds;
d = params.d;

if isfield(params, 'seeds') && isfield(params, 'kdtree')
    kdtree = params.kdtree;
else
    kdtree = vl_kdtreebuild(seeds');
end

if ~isfield(params, 'nn')
    nn = 10;
else
    nn = params.nn;
end

if ~isfield(params, 'useFilter')
    useFilter = 0;
else
    useFilter = params.useFilter;
end

if useFilter == 1
    kdtree_single = vl_kdtreebuild(p');
    [index_single, distance_single] = vl_kdtreequery(kdtree_single, p', seeds', 'NumNeighbors', 1, 'MaxComparisons', 10);
    index_single2 = index_single .* uint32(distance_single < d/2);
    nonzero_single_idx = find(index_single2(:) > 0);
    seeds_filtered = seeds(nonzero_single_idx,:);
else
    nonzero_single_idx = 1:length(seeds);
    seeds_filtered = seeds;
end

kdtree = vl_kdtreebuild(seeds_filtered');
[index, distance] = vl_kdtreequery(kdtree, seeds_filtered', p',  'NumNeighbors', nn, 'MaxComparisons', 100);
index2 = index .* uint32(distance < d/2);

voxelized = zeros(size(seeds, 1)*2, 3);

%For color
c_seed_all = zeros(size(seeds_filtered, 1), 3);
c_seed_cnt = zeros(size(seeds_filtered, 1), 1);
for ii=1:size(index2,1)
    %ii
    nonzero_idx = find(index2(ii,:) > 0);
    c_seeds1 = accumarray(index2(ii,nonzero_idx)', c(nonzero_idx, 1));
    c_seeds2 = accumarray(index2(ii,nonzero_idx)', c(nonzero_idx, 2));
    c_seeds3 = accumarray(index2(ii,nonzero_idx)', c(nonzero_idx, 3));
    
    c_seed_all(1:length(c_seeds1), 1) = c_seed_all(1:length(c_seeds1), 1) + c_seeds1;
    c_seed_all(1:length(c_seeds2), 2) = c_seed_all(1:length(c_seeds2), 2) + c_seeds2;
    c_seed_all(1:length(c_seeds3), 3) = c_seed_all(1:length(c_seeds3), 3) + c_seeds3;
    
    c_seed_cnt(1:length(c_seeds1)) = c_seed_cnt(1:length(c_seeds1), 1)+(c_seeds1>0);
end
nz_idx = find(c_seed_cnt > 0);
c_seed_all(nz_idx, :) = c_seed_all(nz_idx, :) ./ repmat(c_seed_cnt(nz_idx), 1, 3);
voxelized(nonzero_single_idx,:) = c_seed_all;

%For surface normal
if ~isempty(n)
    n(find(isnan(n))) = 0;
    n_seed_all = zeros(size(seeds_filtered, 1), 3);
    n_seed_cnt = zeros(size(seeds_filtered, 1), 1);
    for ii=1:size(index2,1)
        %ii
        nonzero_idx = find(index2(ii,:) > 0);
        n_seeds1 = accumarray(index2(ii,nonzero_idx)', n(nonzero_idx, 1));
        n_seeds2 = accumarray(index2(ii,nonzero_idx)', n(nonzero_idx, 2));
        n_seeds3 = accumarray(index2(ii,nonzero_idx)', n(nonzero_idx, 3));

        n_seed_all(1:length(n_seeds1), 1) = n_seed_all(1:length(n_seeds1), 1) + n_seeds1;
        n_seed_all(1:length(n_seeds2), 2) = n_seed_all(1:length(n_seeds2), 2) + n_seeds2;
        n_seed_all(1:length(n_seeds3), 3) = n_seed_all(1:length(n_seeds3), 3) + n_seeds3;

        n_seed_cnt(1:length(n_seeds1)) = n_seed_cnt(1:length(n_seeds1), 1)+(n_seeds1>0);
    end
    nz_idx = find(n_seed_cnt > 0);
    n_seed_all(nz_idx, :) = n_seed_all(nz_idx, :) ./ repmat(n_seed_cnt(nz_idx), 1, 3);
    voxelized(size(seeds, 1)+nonzero_single_idx,:) = n_seed_all;
end

end

