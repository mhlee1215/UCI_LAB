function [ projectedPoints ] = projectImg2World( cameraMat, points, depth )
%IMG2WORLD Summary of this function goes here
%   Detailed explanation goes here

% if size(points,1) ~= 3 && size(points,2) == 3 && size(points,2) ~= 1
%     points = points';
% end

if size(depth,1) > size(depth,2)
    depth = depth';
end

if size(points,1) == 2
    points = [points; ones(1,size(points,2))];
end

%points = points';

fmat = [ cameraMat.f 0 cameraMat.imSize(2)/2 ; 0 cameraMat.f cameraMat.imSize(1)/2 ; 0 0 1];

%  P = R * X + t       (conversion from world to camera coordinates)
%  p = -P / P.z        (perspective division)
%  p' = f * r(p) * p   (conversion to pixel coordinates)


p = inv(fmat)*double(points);
%alpha = sqrt((depth.^2) ./ (sum(p'.^2, 2) ));
%p([2 3],:) = -p([2 3],:);
%p = -p;
P = (p .* repmat(depth./p(3,:), 3, 1));
projectedPoints = cameraMat.R'*(P - repmat(cameraMat.t, 1, size(P, 2)));
projectedPoints = projectedPoints';




% fmat = [ cameraMat.f 0 0 ; 0 cameraMat.f 0 ; 0 0 1];
% 
% %camMat = [cameraMat.R, cameraMat.t'];
% %camMat = [camMat; 0 0 0 1];
% cameraPos = inv(fmat)*points';
% alpha = sqrt((depth.^2) ./ (sum(cameraPos'.^2, 2) ));
% cameraPos = cameraPos.* repmat(alpha', 3, 1);
% 
% projectedPoints = cameraMat.R'*(cameraPos - repmat(cameraMat.t', 1, size(cameraPos, 2)));
% projectedPoints = projectedPoints';



end
% 

% function [ projectedPoints ] = projectImg2World( cameraMat, points, depth )
% %IMG2WORLD Summary of this function goes here
% %   Detailed explanation goes here
% 
% if length(points) == 2
%     points = [points 1];
% end
% 
% fmat = [ cameraMat.f 0 0 ; 0 cameraMat.f 0 ; 0 0 1];
% 
% cameraPos = inv(fmat)*points';
% alpha = sqrt(depth*depth / (cameraPos' * cameraPos));
% cameraPos = cameraPos.* alpha;
% 
% projectedPoints = cameraMat.R'*(cameraPos - cameraMat.t');
% projectedPoints = projectedPoints';
% 
% 
% end

