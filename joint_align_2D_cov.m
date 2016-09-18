addpath(genpath('libs'));

X = rand(100, 2);
X = [X zeros(size(X, 1), 1)];
Y = rand(100, 2);
Y = [zeros(size(Y, 1), 1) Y];

cSize = 3;
XX = [];
for i = 1:cSize
    XX = [XX ; bsxfun(@plus, X, [0 0 i*2])];
    XX = [XX ; bsxfun(@plus, Y, [i*2 0 0])];
end



% figure;scatter3(XX(:,1), XX(:, 2), XX(:, 3), 8, 'k', 'filled');
YY = XX;
randIdx = randperm(size(XX, 1));
Xin = XX(randIdx(1:cSize*2), :)';
V = {};
V{1} = XX';
V{2} = YY';
V = V';


% dataRoot = '/home/mhlee/data_from_odroid/complete';
% data_set = {};
% % data_set{1} = 'LAB_1-2016-07-18_13_14.klg_cvt.ply';
% % data_set{2} = 'LAB_1-2016-07-29_12_28.klg_cvt.ply';
% data_set{1} = 'LAB_2-2016-07-18_13_23.klg_cvt.ply';
% data_set{2} = 'LAB_2-2016-07-18_13_25.klg_cvt.ply';
% [dataSet, poseSet] = loadEFDataset( dataRoot, data_set);
% dataSetSmall = loadUniformSampling(dataSet, 15);
% 
% colors = distinguishable_colors(30, [0 0 0]);
% idx1 = 1;
% idx2 = 2;
% 
% % load('/home/mhlee/data_from_odroid/cache/LAB_1-2016-07-18_13_14.mat');
% featureNew = dataSetSmall{idx1};
% % validIdx = 1:length(featureNew.v(:,1));%find((featureNew.v(:,1) > 0));% .* (featureNew.v(:,2) > 0) .* (featureNew.v(:,3) > 0));
% validIdx1 = find((featureNew.v(2,:) > -3).*(featureNew.v(2,:) < -1));% .* (featureNew.v(1,:) < 3));% .* (featureNew.v(2,:) > 1));% .* (featureNew.v(2,:) < 2.5));% .* (featureNew.v(1,:) > 0));
% iValidIdx1 = zeros(length(featureNew.v(1,:)), 1);
% iValidIdx1(validIdx1) = 1:length(validIdx1);
% % v = featureNew.v';%(validIdx, :);
% % c = featureNew.c';%(validIdx, :);
% Y = featureNew.v(:, validIdx1)';
% Yc = featureNew.c(:, validIdx1)';
% Yn = featureNew.n(:, validIdx1)';
% Yi = featureNew.sIndex(find(ismember(featureNew.sIndex, validIdx1)));
% YB = dataSet{idx1};
% bigSubIdx1 = find(ismember(featureNew.sIndex, validIdx1));
% YB.vs = YB.v(:, bigSubIdx1);
% YB.cs = YB.c(:, bigSubIdx1);
% YB.ns = YB.n(:, bigSubIdx1);
% 
% % 
% % manupulateIdx1 = find((Y(:,1) < 1));% .* (config.scene(1,:) > 0) .* (config.scene(1,:) > 0));
% % % config.model2 = config.model;
% % Y(manupulateIdx1, 3) = Y(manupulateIdx1, 3)+0.3;% + 0.1.*rand(size(v))-0.1;
% % config.model = config.model2;
% 
% % pclviewer([config.model2 repmat([1 0 0], size(config.scene, 1), 1) ; config.model repmat([1 1 0], size(config.model, 1), 1)]');
% 
% % load('/home/mhlee/data_from_odroid/cache/LAB_1-2016-07-29_12_28.mat');
% featureNew = dataSetSmall{idx2};
% 
% validIdx2 = find((featureNew.v(2,:) > -3).*(featureNew.v(2,:) < -1));% .* (featureNew.v(1,:) < 3));
% iValidIdx2 = zeros(length(featureNew.v(1,:)), 1);
% iValidIdx2(validIdx2) = 1:length(validIdx2);
% X = featureNew.v(:, validIdx2)';% + 0.1.*rand(size(v))-0.1;
% Xc = featureNew.c(:, validIdx2)';% + 0.1.*rand(size(v))-0.1;
% Xn = featureNew.n(:, validIdx2)';% + 0.1.*rand(size(v))-0.1;
% XB = dataSet{idx2};
% Xi = featureNew.sIndex(find(ismember(featureNew.sIndex, validIdx2)));
% bigSubIdx2 = find(ismember(featureNew.sIndex, validIdx2));
% XB.vs = XB.v(:, bigSubIdx2);
% XB.cs = XB.c(:, bigSubIdx2);
% XB.ns = XB.n(:, bigSubIdx2);
% 
% 
% M=2;
% va = [];
% for i = 1:M
%     va = [va V{i}];
% end
% clusterDensity = 20;
% [Xin,~,~,~,density]= uniformSubSample( va', clusterDensity);
% Xin = Xin(:, find(density > 1));
% K = size(Xin, 2);
% 
% 
% V = {Y' ; X'};
% N = {Yn' ; Xn'};
% C = {Yc' ; Xc'};

% [X1 Y1] = meshgrid(1:10, 1:10);
% XX = [X1(:) Y1(:) zeros(length(X1(:)), 1)]';
% V = {};
% V{1} = XX;
% V{2} = XX;
% V = V';
% Xin = mean(XX, 2);

maxNumIter = 150;
updateVis = 0 ;

[R,t,X2,S,a,pk,T,TAssigned, TXQ, vis, ~, ~] = jrmpc_soft_with_cov2d(V,Xin, ...
    'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0, 'updateVis', updateVis, 'updateTR', 0);
% [R,t,X2,S,a,pk,T] = jrmpc(V, Xin,'maxNumIter',maxNumIter,'gamma', 0.1);%, 'updateTR', 0);

figure; scatter3(V{1}(1,:), V{1}(2,:), V{1}(3,:));
hold on;
for i=1:size(X2, 2)
    h1 = plot_gaussian_ellipsoid(X2(:,i), S(:,:,i), 1);
%     h1 = plot_gaussian_ellipsoid(X2(:,i), eye(3).*S(i), 10);
    set(h1,'facealpha',0.2);
end
axis equal;

params.type = 2;
params.Assigned = TAssigned;
params.K = length(TXQ{1, 1, end});
params.interval = 5;
% params.marker = marker;
% params.markerSize = markerSize;
% params.clrmap = clrmap;
% params.strIdx = strIdx;
params.view = [-12 61];
params.pause = 0.1;
params.TXQ = TXQ;
h = figure;
[TV] = drawTransformation(V, T, params);













Vertices = bsxfun(@times, rand(3, 100), [1 3 1]');
X = mean(Vertices,2);
Q = cov(Vertices');
figure; scatter3(Vertices(1,:), Vertices(2,:), Vertices(3,:));
hold on;
h1 = plot_gaussian_ellipsoid(X, Q, mean(mean(Q))*10 );
set(h1,'facealpha',0.2);
axis equal;

Q2 = covDimReduction(Q);
figure; scatter3(Vertices(1,:), Vertices(2,:), Vertices(3,:));
hold on;
h1 = plot_gaussian_ellipsoid(X, Q2, mean(mean(Q))*10 );
set(h1,'facealpha',0.2);
axis equal;
