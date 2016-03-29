function [ validSeeds ] = getValidSeeds( p, seeds, d )
%GETVALIDSEEDS Summary of this function goes here
%   Detailed explanation goes here

    kdtree_single = vl_kdtreebuild(p');
    [index_single, distance_single] = vl_kdtreequery(kdtree_single, p', seeds', 'NumNeighbors', 1, 'MaxComparisons', 10);
    index_single2 = index_single .* uint32(distance_single < d/2);
    nonzero_single_idx = find(index_single2(:) > 0);
    validSeeds = seeds(nonzero_single_idx,:);
    
end

