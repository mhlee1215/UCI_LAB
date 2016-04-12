% t = -2:2;
% dt = 1;
% ti =-2:0.025:2;
% dti = 0.025;
% y = sign (t);
% ys = interp1 (t,y,ti,'spline');
% yp = interp1 (t,y,ti,'pchip');
% ddys = diff (diff (ys)./dti) ./ dti;
% ddyp = diff (diff (yp)./dti) ./ dti;
% figure (1);
% plot (t, y, 'b-', ti,ys,'r-', ti,yp,'g-');
% legend ('spline', 'pchip', 4);
% figure (2);
% plot (ti,ddys,'r+', ti,ddyp,'g*');
% legend ('spline', 'pchip');




p = [[0 0 0];
    [0 0 1];
    [0 1 1];
    [1 1 1];
    [0 0 0]
    ];

dt = 0.01;
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
scatter3(x, y, z, 50, 'r', 'filled');
scatter3(pp(:,1), pp(:,2), pp(:,3), 10, 'b', 'filled');



p1 = pp(1,:);
p2 = pp(2,:);
r = vrrotvec([0 0 1], p2 - p1);
m = vrrotvec2mat(r);



[pp tt] = fnplt(cscvn(points), 200); 
hold on, 
pp = pp';
scatter3(pp(:,1), pp(:,2), pp(:,3), 'filled');
plot(points(1,:),points(2,:),'o'), hold off

p = [[0 0 0];
    [0 0 1];
    [0 1 1];
    [1 1 1];
    [0 0 0]
    ];
curve = cscvn(p);
figure; hold on;
scatter3(p(:,1), p(:,2), p(:,3), 'filled');
fnplt(curve);




x = p(:,1);%0:0.1:1;
y = p(:,2);%x;!
z = p(:,3);%x;
f = @(x,y,z) x.^2 - y - z.^2;
[xx, yy, zz] = meshgrid (x, y, z);
v = f (xx,yy,zz);
xi = 0:0.05:1;
yi = xi;
zi = xi;
[xxi, yyi, zzi] = meshgrid (xi, yi, zi);
vi = interp3 (x, y, z, v, xxi, yyi, zzi, 'spline');
vi = interp3 (p, 'spline');
% [xxi, yyi, zzi] = ndgrid (xi, yi, zi);
% vi2 = interpn (x, y, z, v, xxi, yyi, zzi, 'spline');

mesh (zi, yi, squeeze (vi(1,:,:)));







p1 = [-0.5 1 -3];
p2 = [-0.5 1 -1];


r = vrrotvec([0 0 1], p2-p1);
m = vrrotvec2mat(r);