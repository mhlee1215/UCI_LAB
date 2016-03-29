function [ vertices_sub, index, color_sub ] = uniformSubSample( vertices, qFactor, colors)
%SUBSAMPLE Summary of this function goes here
%   Detailed explanation goes here
% 
%     if isempty(minP)
    color_sub = [];

     minP = min(vertices);
%     end
%     maxP = max(vertices);
% qFactor = 100;
    qVertices = uint16((vertices - repmat(minP, size(vertices, 1), 1)).*qFactor)+1;
    
    sXMean = accumarray(qVertices, vertices(:,1), [], @mean);
    sYMean = accumarray(qVertices, vertices(:,2), [], @mean);
    sZMean = accumarray(qVertices, vertices(:,3), [], @mean);
    
    validIdx = find((sXMean ~= 0).*(sYMean ~= 0).*(sZMean ~= 0));
    vertices_sub = [sXMean(validIdx), sYMean(validIdx), sZMean(validIdx)]';
    
    kdtree = vl_kdtreebuild(vertices_sub) ;
    [index, distance] = vl_kdtreequery(kdtree, vertices_sub, vertices', 'MaxComparisons', 55) ;
    
    if exist('colors', 'var')
        sRMean = accumarray(qVertices, colors(:,1), [], @mean);
        sGMean = accumarray(qVertices, colors(:,2), [], @mean);
        sBMean = accumarray(qVertices, colors(:,3), [], @mean);
    
        color_sub = [sRMean(validIdx), sGMean(validIdx), sBMean(validIdx)]';
    end
    
   
end

