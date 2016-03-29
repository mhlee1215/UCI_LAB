function [ imgPoint ] = projectWorld2Img( cameraMat, points)
%PROJECT Summary of this function goes here
%   Detailed explanation goes here
%   points = N x 3

fmat = [ cameraMat.f 0 cameraMat.imSize(2)/2 ; 0 cameraMat.f cameraMat.imSize(1)/2 ; 0 0 1];
camMat = [cameraMat.R cameraMat.t'];
points = [points ones(size(points, 1), 1)];

P = camMat*points';
%P([2 3],:) = -P([2 3],:);
p = (P ./ repmat(P(3,:), 3, 1));
imgPoint = fmat * p;
imgPoint = imgPoint';
imgPoint = imgPoint(:, 1:2);


% fmat = [ cameraMat.f 0 0 ; 0 cameraMat.f 0 ; 0 0 1];
% %cMat = fmat * camMat;
% 
% points = [points ones(size(points, 1), 1)];
% % camMat = [cameraMat.R cameraMat.t'];
% projectedPoints = fmat*camMat*points';
% %projectedPoints = cMat*points';
% projectedPoints = projectedPoints';
% 
% projectedPoints(:,1) = projectedPoints(:,1) ./ projectedPoints(:,3);
% projectedPoints(:,2) = projectedPoints(:,2) ./ projectedPoints(:,3);
% 
% projectedPoints(:,3) = projectedPoints(:,3) ./ projectedPoints(:,3);
% 
% imgPoint = projectedPoints(:, 1:2);

end

% function [ imgPoint ] = projectWorld2Img( cameraMat, points )
% %PROJECT Summary of this function goes here
% %   Detailed explanation goes here
% 
% 
% fmat = [ cameraMat.f 0 0 ; 0 cameraMat.f 0 ; 0 0 1];
% %cMat = fmat * camMat;
% 
% points = [points ones(size(points, 1), 1)];
% camMat = [cameraMat.R cameraMat.t'];
% projectedPoints = fmat*camMat*points';
% %projectedPoints = cMat*points';
% projectedPoints = projectedPoints';
% 
% projectedPoints(:,1) = projectedPoints(:,1) ./ projectedPoints(:,3);
% projectedPoints(:,2) = projectedPoints(:,2) ./ projectedPoints(:,3);
% 
% %projectedPoints(:,3) = projectedPoints(:,3) ./ projectedPoints(:,3);
% 
% imgPoint = projectedPoints(:, 1:2);
% 
% end