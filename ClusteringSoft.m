function [ groups ] = softClustering( fv1_v, gKey1 )
%SOFTCLUSTERING Summary of this function goes here
%   Detailed explanation goes here

    kdtree = vl_kdtreebuild(gKey1) ;
    [index, distance] = vl_kdtreequery(kdtree, gKey1, fv1_v', 'NumNeighbors', 5, 'MaxComparisons', 55) ;

    index_filtered = index;
    %Filter out by distance
    filterout_idx = find(distance >= max(distance(1, :)));
    index_filtered(filterout_idx) = 0;
    %Filter out by rank
    topN = 3;
    rankMat = zeros(size(distance));
    rankMat(topN+1:end, :) = 1;
    filterout_idx = find(rankMat);
    index_filtered(filterout_idx) = 0;

    groups = {};
    for i=1:size(gKey1, 2)
        groups{i} = find(sum(index_filtered == i, 1));
    end

end

