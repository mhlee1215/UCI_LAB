function [ index_filtered, index, distance ] = getKNN( C_src, C_dst )
%GETKNN Summary of this function goes here
%   Detailed explanation goes here

% C_src = clustersCenters{g_id_src};
% C_dst = clustersCenters{g_id_dst};


kdtree = vl_kdtreebuild(C_dst) ;
[index, distance] = vl_kdtreequery(kdtree, C_dst, C_src, 'NumNeighbors', 10, 'MaxComparisons', 55) ;

index_filtered = index;
filterout_idx = (distance >= max(distance(1, :)*0.95));
index_filtered(filterout_idx) = 0;

end

