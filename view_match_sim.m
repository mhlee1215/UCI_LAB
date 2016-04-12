function [ output_args ] = view_match_sim( idx, interval, DescrData, vertex11, color1, vertex22, color2, key11, key22, NN_Data, searchRadius )
%VIEW_MATCH_SIM Summary of this function goes here
%   Detailed explanation goes here

topK = 5;

asub = DescrData(idx,:);
validIdx = find(asub < 500);
invalidIdx = find(asub >= 500);
if ~isempty(validIdx)
    asub = asub - min(asub(validIdx));
    asub = asub ./ max(asub(validIdx));
end
asub(invalidIdx) = 1;
color_zeros = ones(length(asub), 1);
color = [color_zeros (1-asub)' (1-asub)'];

scatter3(vertex11(1:interval:end,1), vertex11(1:interval:end,2), ...
    vertex11(1:interval:end,3), 8, color1(1:interval:end, :), 'filled');

scatter3(vertex22(1:interval:end,1), vertex22(1:interval:end,2), ...
    vertex22(1:interval:end,3), 8, color2(1:interval:end, :), 'filled');

scatter3(key11(1,idx)', key11(2,idx)', key11(3,idx)', 200, 'filled');
drawSphere(key11(1,idx)', key11(2,idx)', key11(3,idx)', searchRadius);
scatter3(key22(1,:)', key22(2,:)', key22(3,:)', 80, color, 'filled');

% scatter3(key22(1,NN_Data(idx, 1:topK))', key22(2,NN_Data(idx, 1:topK))', key22(3,NN_Data(idx, 1:topK))', 160, 'k', 'filled');



end

