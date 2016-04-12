function [ poseSet ] = syn_path_gen( p, dt, vertices )
%SYN_PATH_GEN Summary of this function goes here
%  in : N x 3 Mat


% dt = 0.01;
curve = p';
% interpote a 3d curve using spline 
% path 3*x 
% newPath 3*x

x = curve(1, :); 
y = curve(2, :); 
z = curve(3, :);

t = cumsum([0;sqrt(diff(x(:)).^2 + diff(y(:)).^2 + diff(z(:)).^2)]); 
sx = spline(t,x); 
sy = spline(t,y); 
sz = spline(t,z);

tt = t(1):dt:t(end); 
xp = ppval(sx, tt); 
yp = ppval(sy, tt); 
zp = ppval(sz, tt);

newCurve = [xp; yp; zp]; 
pp = newCurve';
figure; hold on;
scatter3(x, y, z, 100, 'o', 'filled');

scatter3(x(1), y(1), z(1), 50, 'r', 'filled');
scatter3(x(end), y(end), z(end), 50, 'g', 'filled');

scatter3(pp(:,1), pp(:,2), pp(:,3), 10, 'b', 'filled');


scatter3(vertices(1:10:end,1),vertices(1:10:end,2), vertices(1:10:end,3), 0.5, 'filled');

poseSet = {};
for i=1:size(pp, 1)-1

    p1 = pp(i,:);
    p2 = pp(i+1,:);
    r = vrrotvec([0 0 -1], -(p2 - p1));
    m = vrrotvec2mat(r);
    
%     m(2,2) = m(2,2)*-1;
    m = m.*[1 1 1 ; -1 -1 -1; 1 1 1];
    
    zz = eye(3);
    
%     zz(1,1) = zz(1,1)*-1;
%     zz(3,3) = zz(3,3)*-1;
    
    m = inv(m)*zz;
    
    poseSet{end+1}.R = m;
    poseSet{end}.t = -m*p1';
    
end

% imW = 1;
% imH = 1;
% for i=1:1:length(poseSet)
%     cameraMat = poseSet{i};
% %    pM(:,end)
% %     scatter3(cameraMat.t(1), cameraMat.t(2), cameraMat.t(3), 0.1);
% %     P = projectImg2World(pM, [0 0]', 0.1);
%     
%     cameraMat.f = 1;
%     cameraMat.imSize = [0 0];%[imH imW];
% %     cameraMat.R
%     Cp = -cameraMat.R'*cameraMat.t;
%    
%     depth = 0.5;
%     C1 = projectImg2World(cameraMat, [-imH -imW]', depth);
%     C2 = projectImg2World(cameraMat, [-imH imW]', depth);
%     C3 = projectImg2World(cameraMat, [imH imW]', depth);
%     C4 = projectImg2World(cameraMat, [imH -imW]', depth);
%     PScene = [C1 ; C2 ; C3 ; C4 ; C1];
% 
% %     scatter3(Cp(1), Cp(2), Cp(3), 50, 'm', 'filled');
%     for i=1:size(PScene, 1)-1
%         line([Cp(1) PScene(i,1)]', [Cp(2) PScene(i,2)]', [Cp(3) PScene(i,3)]', 'linewidth', 1);    
%         line([PScene(i,1) PScene(i+1,1)]', [PScene(i,2) PScene(i+1,2)]', [PScene(i,3) PScene(i+1,3)]', 'linewidth', 1);    
%     end
%     
% end
axis equal;
grid;





end

