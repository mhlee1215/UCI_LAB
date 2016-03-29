function [ T ] = GMBasedRegistration( v_m1_s, c_m1_s, v_m2_s, c_m2_s, colorLambda )
%GMBASEDREGISTRATION Summary of this function goes here
%   Detailed explanation goes here

V = {};
V{end+1} = v_m1_s';
V{end+1} = v_m2_s';
V2 = V';

C = {};
C{end+1} = c_m1_s';
C{end+1} = c_m2_s';
C2 = C';

d_size = 2;
maxNumIter = 100;            
M = d_size;
idx = transpose(1:M);
strIdx = arrayfun(@(j) sprintf('view %d',j),idx,'uniformoutput',false);

K = ceil(0.5*median(cellfun(@(V) size(V,2),V))); 
K = 1000;

% rGMM_color = bsxfun(@plus, rand(K, 3)./2, [0.3 0.3 0.3]);

% sample the unit sphere, by randomly selecting azimuth / elevation angles
% az = 2*pi*rand(1,K);
% el = 2*pi*rand(1,K);
% 
% %points on a unit sphere
% Xin = [cos(az).*cos(el); sin(el); sin(az).*cos(el)];% (unit) polar to cartesian conversion
% 
% Xin = Xin/10; % it is good for the initialization to have initial cluster centers at the same order with the points
% % since sigma is automatically initialized based on X and V
% 
% meanData = cellfun(@(a) mean(a')',V,'uniformoutput',false); % 1 x K rows
% meanData = mean(cat(3,meanData{:}),3);
% % Xin = Xin ./ repmat(sqrt(var(Xin'))', 1, size(Xin, 2));
% Xin = Xin + repmat(meanData, 1, size(Xin, 2));

VV = V2{1};
Xin = zeros(3, K);
randIdx = randperm(size(VV, 2));
Xin(:,:) = VV(:, randIdx(1:K));

if colorLambda == 0.0
    [R,t,X,S,a,pk,T,TAssigned, TXQ] = jrmpc_soft(V2,Xin, ...
    'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0);
else
    [R,t,X,S,a,pk,T,TAssigned, TXQ] = jrmpc_soft_color(V2,Xin, ...
    'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0, ...
    'color', C2, 'colorLambda', colorLambda);
end





markerSize = {};
for ii=1:d_size
    markerSize{end+1} = 1+ii;
end
markerSize = markerSize';
markerSize = markerSize(1:d_size);

% markers
marker_set = {'s', 'x', '.', '^'};
marker = {};
for i=1:d_size
    marker{end+1} = marker_set{mod(i, 4)};
end
marker = marker';

clrmap = {[1 1 0]; [1 0 1] ; [0 1 1] ; [1 0 0] ; [0 1 0] ; [0 0 1] ; [1 1 1] ; [0 0 0] ; [0.3 0.3 0.3] ; [.5 0.9 0.9] ; [0.8 1 1]; [0 0 0] ...
    ; [1 1 0]; [1 0 1] ; [0 1 1] ; [1 0 0] ; [0 1 0] ; [0 0 1] ; [1 1 1] ; [0 0 0] ; [0.3 0.3 0.3] ; [.5 0.9 0.9] ; [0.8 1 1]; [0 0 0]};
clrmap = clrmap(1:d_size);

params.type = 2;
params.Assigned = TAssigned;
params.K = K;
params.interval = 3;
params.marker = marker;
params.markerSize = markerSize;
params.clrmap = clrmap;
params.strIdx = strIdx;
params.view = [];%[40 54];
params.pause = .1;
h = figure;
[TV] = drawTransformation(V2, T, params);
axis equal;

end

