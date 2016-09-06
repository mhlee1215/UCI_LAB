
% [xx,yy] = meshgrid(-10:0.2:10,-10:0.2:10);
% %zz = 5*sin(yy/5)+0.1*(xx.^2+yy.^2).^0.5;
% zz = 0.5*(10-sqrt(5+xx.^2+yy.^2))+ 2*sin(xx/5);
% zz = zz + 0.1*randn(size(zz));
% zz = zz - min(zz(:));

[xx yy zz] = sphere(80);

p1 = [xx(:) yy(:) zz(:)];
c1 = ones(3, length(xx(:)));

idx = find((p1(:,1) > 0) .* (p1(:,2) > 0) .* (p1(:,3) > 0));
c1(1,idx) = 0;
idx2 = find((p1(:,1) < 0) .* (p1(:,2) < 0) .* (p1(:,3) < 0));
c1(2:3,idx2) = 0;

noise = 0.03;
R = eye(3);
z = pi/4;
R(1,1) = cos(z);
R(1,2) = -sin(z);
R(2,1) = sin(z);
R(2,2) = cos(z);

idx1 = find(p1(:,1) < 0.4);
idx2 = find(p1(:,1) > -0.4);


p2 = p1(idx2, :)';
c2 = c1(:,idx2);%ones(3, length(xx(:)));
p22 = bsxfun(@plus, rand(3, 20000).*0.7, [2 2 0]');
c22 = ones(size(p22));


p2 = bsxfun(@plus, R*p2 , [0.3 0.3 0 ]')+noise*randn(size(p2));
p2 = [p2 p22];
c2 = [c2 c22];

p1 = p1(idx1, :);
p1 = p1 + noise*randn(size(p1));

c1 = c1(:, idx1);
p11 = bsxfun(@plus, rand(3, 20000).*0.7, [2 2 0]');
c11 = ones(size(p22));
p1 = [p1 ;p11'];
c1 = [c1' ; c11']';


% pclviewer([p1 c1' ; p2' c2']');

V2 = {};
V2{1} = p1';
V2{2} = p2;
C2 = {};
C2{1} = c1;
C2{2} = c2;
V2=V2';
C2=C2';

[R2,t2,X2,S2,a2,pk2,T2,TAssigned2, TXQ2, vis2] = joint_align(V2,[], C2, 100);

% 
TV2 = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V2,T2(:,1,end),T2(:,2,end),'uniformoutput',false);
%                 TVb2 = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),Vb,T2(:,1,end),T2(:,2,end),'uniformoutput',false);
pclviewer([TV2{1} TV2{2}; c1 c2]);

params.type = 1;
params.Assigned = TAssigned2;
params.K = length(TXQ2{1, 1, end});
params.interval = 2;
params.view = [];%[40 54];
params.pause = 0.1;
params.TXQ = TXQ2;
h = figure;
drawTransformation(V2, T2, params);