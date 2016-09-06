function [ pointsNorm ] = fcn_floorNormalizationCharless( points, pointsSeed, inlierIndex )
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

xx = pointsSeed(inlierIndex, 1)';
yy = pointsSeed(inlierIndex, 2)';
zz = pointsSeed(inlierIndex, 3)';
st = tpaps([xx;yy],zz);


% figure(1); clf;
% fnplt(st), hold on
% % plot3(xx(p),yy(p),zz(p),'wo','markerfacecolor','k')
% plot3(xx,yy,zz,'wo','markerfacecolor','k')
% axis image;
% axis vis3d;

gsize = 0.5;
[xx2,yy2] = meshgrid(min(xx):gsize:max(xx),min(yy):gsize:max(yy));
zz2 = fnval(st,[xx2(:) yy2(:)]');
zz2 = reshape(zz2,size(xx2,1),size(xx2,2));

% estimate surface normal using interpolated grid
[dxi,dyi] = gradient(zz2,gsize,gsize);
dx = -dxi ./ sqrt(dxi.^2 + dyi.^2 + 1);
dy = -dyi ./ sqrt(dxi.^2 + dyi.^2 + 1);
dz = ones(size(dx)) ./ sqrt(dxi.^2 + dyi.^2 + 1);


% original points
src1 = [xx2(:)'; yy2(:)'; zz2(:)']';
target1 = [xx2(:)'; yy2(:)'; zeros(size(zz2(:)'))]';

% points offset along surface normal
sc = max(pointsSeed(:, 3)); %size of offset
src2 = [xx2(:)'+sc*dx(:)'; yy2(:)'+sc*dy(:)'; zz2(:)' + sc*dz(:)']';
target2 = [xx2(:)'+dx(:)'; yy2(:)'+dy(:)'; sc*dz(:)']';



% plot normals
% figure(1);
% hold on;
% for i = 1:size(src1, 1);
%   plot3([src1(i,1) src2(i,1)],[src1(i,2) src2(i,2)],[src1(i,3) src2(i,3)],'-');
% end



% src = pointsSeed(inlierIndex(1:10:end), :);
% target = src;
% target(:, 3) = 0;
% 
% maxZ = max(pointsSeed(:, 3));
% src = [src ; bsxfun(@plus, src, [0 0 maxZ])];
% target = [target ; bsxfun(@plus, target, [0 0 maxZ])];

src = [src1; src2];
target = [target1; target2];

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

