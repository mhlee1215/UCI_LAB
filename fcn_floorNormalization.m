function [ pointsNorm ] = fcn_floorNormalization( points, pointsSeed, inlierIndex )
%FCN_FLOORNORMALIZATION Summary of this function goes here
%   Detailed explanation goes here
% v = (R*points')';
% zOffset = mean(Rv(inliers, 3));
% T = -[0 0 zOffset]';
% 
% pclviewer([bsxfun(@plus, Rv',T)]);

% points2 = bsxfun(@plus, Rv',T)';
% points2 = bsxfun(@plus, Rv',T)';

% warped = fcn_floorNormalization(points, inliers);

src = pointsSeed(inlierIndex(1:10:end), :);
target = src;
target(:, 3) = 0;

maxZ = max(pointsSeed(:, 3));
src = [src ; bsxfun(@plus, src, [0 0 maxZ])];
target = [target ; bsxfun(@plus, target, [0 0 maxZ])];

pointsNorm = zeros(size(points));
step = 100000;
curStep = 1;
    disp(sprintf('Do floor normalization..%d pts', size(points, 1)));
while curStep < size(points, 1)
%     disp(sprintf('%d/%d', curStep, size(points, 1)));
    nextStep = min(curStep + step, size(points, 1));
    pointsNorm(curStep:nextStep, :) = warp3Dtps2(src,target,points(curStep:nextStep, :), 1);    
    curStep=nextStep;
end


% vertices_sub = seedPoints(inlierIndex, :);
% vertices = points;
% kdtree = vl_kdtreebuild(vertices_sub') ;
% [index, distance] = vl_kdtreequery(kdtree, vertices_sub', vertices', 'MaxComparisons', 55) ;
% 
% pointsNorm = points;
% pointsNorm(:, 3) = pointsNorm(:, 3) - seedPoints(inlierIndex(index), 3);

end

