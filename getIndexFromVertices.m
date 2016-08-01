function [ SegmentIndex, segmentId ] = getIndexFromVertices( V, Segments )
%GETINDEXFROMVERTICES Summary of this function goes here
%   Detailed explanation goes here

disp(sprintf('Build kd-tree...\n'));
kdtree = vl_kdtreebuild(V');
SegmentIndex = cell(1, length(Segments));
segmentId = zeros(size(V, 1), 1);
for i=1:length(Segments)
    disp(sprintf('Get indices for segment %d/%d', i, length(Segments)));
    [index, ~] = vl_kdtreequery(kdtree, V', double(Segments{i}(1:3,:)), 'MaxComparisons', 10) ;
    SegmentIndex{i} = index;
    segmentId(index) = i;
end




end

