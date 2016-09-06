function [R, T, inliers] = fcn_convertGroundNormalCore( points, normals, colors )
%FCN_CONVERTGROUNDNORMAL Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3; colors = zeros(size(points));  end;

disp('Ransac for fit plane...');
% progressbar2(curStep/maxStep); curStep = curStep + 1;
[B, ~, inliers] = ransacfitplane(points', 0.07, colors');




% pclviewer([points(inliers,:)]');
baseNormal = B(1:3);%mean(nForSeg(inliers, :));
meanNormal = mean(normals(inliers, :));

angle = acos(baseNormal'*meanNormal');

if angle > pi/2
    baseNormal = -baseNormal;
end
% if baseNormal(1) * meanNormal(1) < 0
%     baseNormal = baseNormal .* -1;
% end

R=fcn_RotationFromTwoVectors(baseNormal, [0 0 1]);


Rv = (R*points')';
zOffset = mean(Rv(inliers, 3));
T = -[0 0 zOffset]';

% pclviewer([bsxfun(@plus, Rv',T)]);
% 
% points2 = bsxfun(@plus, Rv',T)';
% % points2 = bsxfun(@plus, Rv',T)';
% 
% % warped = fcn_floorNormalization(points, inliers);
% vertices_sub = points2(inliers, :);
% vertices = points2;
% kdtree = vl_kdtreebuild(vertices_sub') ;
% [index, distance] = vl_kdtreequery(kdtree, vertices_sub', vertices', 'MaxComparisons', 55) ;
% 
% points3 = points2;
% points3(:, 3) = points3(:, 3) - points3(inliers(index), 3);
% pclviewer([points3]');

% warped = warp3Dtps2(Rv2(inliers,:),Rv2(inliers,:),points,100);
% 
% pclviewer([Rv2T]');
% pclviewer([warped]');
% pclviewer([points]');

end

