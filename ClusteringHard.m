function [ groups ] = ClusteringHard( fv2_v, gKey2 )
%MYCLUSTERING Summary of this function goes here
%   Detailed explanation goes here


% gKey2 = double(key2(:, pos_NN_Data(validIdx, 1)));
% fv2_v = fv{2}.Vertices;
% fv2_v = bsxfun(@plus, bestR*fv2_v', bestT)';
% fv2_c = fv{2}.FaceVertexCData;
IDX_2 = knnsearch(gKey2', fv2_v);

groups = {};
for i=1:size(gKey2, 2)
    groups{i} = find(IDX_2==i)';
end

end

