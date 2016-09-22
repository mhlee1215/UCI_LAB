
% function [] = EF_Nonrigid_Test()
%  Nonrigid Example 6. Coherent Point Drift (CPD).
%  Nonrigid registration of 3D bunny point sets with use of FGT and Lowrank kernel approximation.

clear all; close all; clc;
%load cpd_data3D_bunny.mat % load bunny set (8171x3) res with defformation

addpath(genpath('./libs'));

% config = gmmreg_load_config('./face.ini');

dataRoot = '/home/mhlee/data_from_odroid/complete';
data_set = {};
% data_set{1} = 'LAB_1-2016-07-18_13_14.klg_cvt.ply';
% data_set{2} = 'LAB_1-2016-07-29_12_28.klg_cvt.ply';
data_set{1} = 'LAB_2-2016-07-18_13_23.klg_cvt.ply';
data_set{2} = 'LAB_2-2016-07-18_13_25.klg_cvt.ply';
data_set{3} = 'LAB_2-2016-07-18_13_26.klg_cvt.ply';
[dataSet, poseSet] = loadEFDataset( dataRoot, data_set);
dataSetSmall = loadUniformSampling(dataSet, 20);

colors = distinguishable_colors(30, [0 0 0]);
idx1 = 1;
idx2 = 2;
idx3 = 3;

% load('/home/mhlee/data_from_odroid/cache/LAB_1-2016-07-18_13_14.mat');
featureNew = dataSetSmall{idx1};
% validIdx = 1:length(featureNew.v(:,1));%find((featureNew.v(:,1) > 0));% .* (featureNew.v(:,2) > 0) .* (featureNew.v(:,3) > 0));
validIdx1 = find((featureNew.v(2,:) > -3).*(featureNew.v(2,:) < -1));% .* (featureNew.v(1,:) < 3));% .* (featureNew.v(2,:) > 1));% .* (featureNew.v(2,:) < 2.5));% .* (featureNew.v(1,:) > 0));
iValidIdx1 = zeros(length(featureNew.v(1,:)), 1);
iValidIdx1(validIdx1) = 1:length(validIdx1);
% v = featureNew.v';%(validIdx, :);
% c = featureNew.c';%(validIdx, :);
Y = featureNew.v(:, validIdx1)';
Yc = featureNew.c(:, validIdx1)';
Yn = featureNew.n(:, validIdx1)';
Yi = featureNew.sIndex(find(ismember(featureNew.sIndex, validIdx1)));
YB = dataSet{idx1};
bigSubIdx1 = find(ismember(featureNew.sIndex, validIdx1));
YB.vs = YB.v(:, bigSubIdx1);
YB.cs = YB.c(:, bigSubIdx1);
YB.ns = YB.n(:, bigSubIdx1);

% 
% manupulateIdx1 = find((Y(:,1) < 1));% .* (config.scene(1,:) > 0) .* (config.scene(1,:) > 0));
% % config.model2 = config.model;
% Y(manupulateIdx1, 3) = Y(manupulateIdx1, 3)+0.3;% + 0.1.*rand(size(v))-0.1;
% config.model = config.model2;

% pclviewer([config.model2 repmat([1 0 0], size(config.scene, 1), 1) ; config.model repmat([1 1 0], size(config.model, 1), 1)]');

% load('/home/mhlee/data_from_odroid/cache/LAB_1-2016-07-29_12_28.mat');
featureNew = dataSetSmall{idx2};

validIdx2 = find((featureNew.v(2,:) > -3).*(featureNew.v(2,:) < -1));% .* (featureNew.v(1,:) < 3));
iValidIdx2 = zeros(length(featureNew.v(1,:)), 1);
iValidIdx2(validIdx2) = 1:length(validIdx2);
X = featureNew.v(:, validIdx2)';% + 0.1.*rand(size(v))-0.1;
Xc = featureNew.c(:, validIdx2)';% + 0.1.*rand(size(v))-0.1;
Xn = featureNew.n(:, validIdx2)';% + 0.1.*rand(size(v))-0.1;
XB = dataSet{idx2};
Xi = featureNew.sIndex(find(ismember(featureNew.sIndex, validIdx2)));
bigSubIdx2 = find(ismember(featureNew.sIndex, validIdx2));
XB.vs = XB.v(:, bigSubIdx2);
XB.cs = XB.c(:, bigSubIdx2);
XB.ns = XB.n(:, bigSubIdx2);



featureNew = dataSetSmall{idx3};
% validIdx = 1:length(featureNew.v(:,1));%find((featureNew.v(:,1) > 0));% .* (featureNew.v(:,2) > 0) .* (featureNew.v(:,3) > 0));
validIdx3 = find((featureNew.v(2,:) > -3).*(featureNew.v(2,:) < -1));% .* (featureNew.v(1,:) < 3));% .* (featureNew.v(2,:) > 1));% .* (featureNew.v(2,:) < 2.5));% .* (featureNew.v(1,:) > 0));
iValidIdx3 = zeros(length(featureNew.v(1,:)), 1);
iValidIdx3(validIdx3) = 1:length(validIdx3);
% v = featureNew.v';%(validIdx, :);
% c = featureNew.c';%(validIdx, :);
Z = featureNew.v(:, validIdx3)';
Zc = featureNew.c(:, validIdx3)';
Zn = featureNew.n(:, validIdx3)';
Zi = featureNew.sIndex(find(ismember(featureNew.sIndex, validIdx3)));
ZB = dataSet{idx3};
bigSubIdx3 = find(ismember(featureNew.sIndex, validIdx3));
ZB.vs = ZB.v(:, bigSubIdx3);
ZB.cs = ZB.c(:, bigSubIdx3);
ZB.ns = ZB.n(:, bigSubIdx3);





V = {Y' ; X' ; Z'};
N = {Yn' ; Xn' ; Zn'};
C = {Yc' ; Xc' ; Zc'};

clusterDensity = 10;
va = [];
for i = 1:length(V)
    va = [va V{i}];
end
[Xin,~,~,~,density]= uniformSubSample( va', clusterDensity);
Xin = Xin(:, find(density > 1));
K = size(Xin, 2);

[v_one,~,~,~,density]= uniformSubSample( va', 20);
v_one = v_one(:, find(density > 0));
% K = size(Xin, 2);

%100 20 generates good results
% [R,t,XX,S,a,pk,T,TAssigned, TXQ, vis, XXc, XXn, Xin] = joint_align(V,N,C, 100, 10);

% [~,~,XX,S,a,pk,T,TAssigned, TXQ, vis, ~, ~] = jrmpc_soft_with_normal(V,Xin, ...
%     'maxNumIter',100,'gamma',0.1, 'updatepriors', 0, 'updateVis', 0, 'updateTR', 0);


% [R,t,XX,S,a,pk,T,TAssigned, TXQ, vis, ~, ~] = jrmpc_soft_with_cov2d(V,Xin, ...
%     'maxNumIter',10,'gamma',0.1, 'updatepriors', 0, 'updateVis', 0, 'updateTR', 1);

V_one = {v_one};
[~,~,XX,S2,a2,pk2,T2,TAssigned2, TXQ2, vis2, ~, ~] = jrmpc_soft_with_cov2d(V_one,Xin, ...
    'maxNumIter',30,'gamma',0.001, 'updatepriors', 0, 'updateVis', 0, 'updateTR', 0);

% [~,~,XX3,S3,a3,pk3,T3] = jrmpc(V, Xin,'maxNumIter',100,'gamma', 0.1);%, 'updateTR', 0);

% TV = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V,T(:,1,end),T(:,2,end),'uniformoutput',false);

% sum(sum(mean(abs(S2), 3)))
figure; hold on;
% for i=1:length(TV)
%     scatter3(TV{i}(1,:), TV{i}(2,:), TV{i}(3,:), 8, colors(i,:), 'filled');
% end

scatter3(v_one(1,:), v_one(2,:), v_one(3,:), 8, colors(1,:), 'filled');
% scatter3(XX2(1,:), XX2(2,:), XX2(3,:), 15, 'm', 'filled');
scatter3(XX(1,:), XX(2,:), XX(3,:), 15, 'm', 'filled');
% scatter3(XX3(1,:), XX3(2,:), XX3(3,:), 15, 'm', 'filled');
% scatter3(Xin(1,:), Xin(2,:), Xin(3,:), 15, 'g', 'filled');
axis equal;

for i=1:1:size(XX, 2)
    vw = 3;
    if size(S2, 2) == 1
        h1 = plot_gaussian_ellipsoid(XX(:,i), eye(3).*S2(i), vw);
%         h1 = plot_gaussian_ellipsoid(XX3(:,i), eye(3).*S2(i), vw);
    else
        h1 = plot_gaussian_ellipsoid(real(XX(:,i)), real(S2(:,:,i)), vw);
        
    end
%     h1 = plot_gaussian_ellipsoid(X2(:,i), eye(3).*S(i), 1);
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



% Init full set of options %%%%%%%%%%
% opt.method='nonrigid'; % use nonrigid registration with lowrank kernel approximation
% opt.numeig=30;                 % leave only 30 larges (out of 8171) eigenvectors/values to approximate G
% opt.eigfgt=0;                  % use FGT to find the largest eigenvectore/values 
% 
% opt.beta=1;            % the width of Gaussian kernel (smoothness)
% opt.lambda=500;          % regularization weight
% 
% opt.viz=1;              % show every iteration
% opt.outliers=0.1;       % use 0.7 noise weight
% opt.fgt=0;              % use FGT to compute matrix-vector products (2 means to switch to truncated version at the end, see cpd_register)
% opt.normalize=0;        % normalize to unit variance and zero mean before registering (default)
% opt.corresp=0;          % compute correspondence vector at the end of registration (not being estimated by default)
% 
% opt.max_it=30;         % max number of iterations
% opt.tol=1e-7;           % tolerance
% opt.sigma2 = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Init full set of options %%%%%%%%%%
opt.method='nonrigid'; % use nonrigid registration with lowrank kernel approximation
opt.numeig=30;                 % leave only 30 larges eigenvectors/values to approximate G
opt.eigfgt=0;                  %  do not use FGT to find the largest eigenvectore/values 

opt.beta=2;            % the width of Gaussian kernel (smoothness)
opt.lambda=3;          % regularization weight

opt.viz=1;              % show every iteration
opt.outliers=0.1;       % noise weight
opt.fgt=0;              % do not use FGT to compute matrix-vector products (2 means to switch to truncated version at the end, see cpd_register)
opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
opt.corresp=0;          % compute correspondence vector at the end of registration (not being estimated by default)

opt.max_it=100;         % max number of iterations
opt.tol=1e-5;           % tolerance
opt.sigma2 = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% data.X = XX';
% data.Y = Y;
% data.Xc = [];%Xc;
% data.Yc = [];%Yc;
% data.Xn = [];%Xn;
% data.Yn = [];%Yn;
% data.cLambda = 0;%0.5;
% data.nLambda = 0;%0.5;

% data.X = XX';
% data.Y = Y;
% data.Xc = XXc';
% data.Yc = Yc;
% data.Xn = XXn';
% data.Yn = Yn;
% data.cLambda = 1;
% data.nLambda = 1;
% opt.auxData = data;
opt.auxData = [];

scale = 1/max(max(S2(1,1,:)), max(max(S2(2,2,:)), max(S2(3,3,:))));

opt.S = S2.*scale.*0.2 + repmat(eye(3).*0.8, 1, 1, size(S2, 3));
close all;
tic;
[TransformY, ~]=cpd_register(XX',Y, opt);
toc;

% data.X = XX';
% data.Y = X;
% data.Xc = XXc';
% data.Yc = Xc;
% data.Xn = XXn';
% data.Yn = Xn;
% data.cLambda = 1;
% data.nLambda = 1;
% opt.auxData = data;
opt.auxData = [];
close all;
tic;
[TransformX, ~]=cpd_register(XX',X, opt);
toc;


opt.auxData = [];
close all;
tic;
[TransformZ, ~]=cpd_register(XX',Z, opt);
toc;

% figure,cpd_plot_iter(X, Y); title('Before');
% figure,cpd_plot_iter(X, Transform.Y);  title('After registering Y to X');
% 

fmY = TransformY.Y;
diffY = (fmY - Y)';
fmX = TransformX.Y;
diffX = (fmX - X)';
fmZ = TransformZ.Y;
diffZ = (fmZ - Z)';


% YBmove = zeros(size(YB.vs));
YBmove = diffY(:, iValidIdx1(Yi));
XBmove = diffX(:, iValidIdx2(Xi));
ZBmove = diffZ(:, iValidIdx3(Zi));
pclviewer([YB.vs+YBmove XB.vs+XBmove ZB.vs+ZBmove; repmat([1 0 0]', 1, size(YB.vs, 2)) repmat([1 1 0]', 1, size(XB.vs, 2)) repmat([0.5 1 0.5]', 1, size(ZB.vs, 2))]);
pclviewer([YB.vs+YBmove XB.vs+XBmove ZB.vs+ZBmove; YB.cs XB.cs ZB.cs]);
pclviewer([YB.vs XB.vs ZB.vs; YB.cs XB.cs ZB.cs]);


pclviewer([YB.vs+YBmove XB.vs+XBmove; repmat([1 0 0]', 1, size(YB.vs, 2)) repmat([1 1 0]', 1, size(XB.vs, 2))]);

dataOutRoot = '/media/mhlee/14A2ACBFA2ACA6A8/E/lab/work/2015_4d_seg/codes/libs/CPD2/examples';
fileName = 'X_XX_Y_seg';
fcn_saveUniformSizeModel( [YB.vs+YBmove XB.vs+XBmove]', [repmat([1 0 0]', 1, size(YB.vs, 2)) repmat([1 1 0]', 1, size(XB.vs, 2))]', [], dataOutRoot, fileName, 0 );
fileName = 'X_XX_Y_seg_ori';
fcn_saveUniformSizeModel( [YB.vs XB.vs]', [repmat([1 0 0]', 1, size(YB.vs, 2)) repmat([1 1 0]', 1, size(XB.vs, 2))]', [], dataOutRoot, fileName, 0 );
fileName = 'X_XX_Y_color';
fcn_saveUniformSizeModel( [YB.vs+YBmove XB.vs+XBmove]', [YB.cs XB.cs]', [], dataOutRoot, fileName, 0 );
fileName = 'X_XX_Y_color_ori';
fcn_saveUniformSizeModel( [YB.vs XB.vs]', [YB.cs XB.cs]', [], dataOutRoot, fileName, 0 );

return;

opt.auxData = [];
tic;
[TransformXY, ~]=cpd_register(X,Y, opt);
toc;

fmXY = TransformXY.Y;
diffXY = (fmXY - Y)';
XYBmove = diffXY(:, iValidIdx1(Yi));
pclviewer([YB.vs+XYBmove XB.vs; repmat([1 0 0]', 1, size(YB.vs, 2)) repmat([1 1 0]', 1, size(XB.vs, 2))]);
pclviewer([YB.vs+XYBmove XB.vs; YB.cs XB.cs]);
pclviewer([YB.vs XB.vs; YB.cs XB.cs]);

dataOutRoot = '/media/mhlee/14A2ACBFA2ACA6A8/E/lab/work/2015_4d_seg/codes/libs/CPD2/examples';
fileName = 'XY_seg';
fcn_saveUniformSizeModel( [YB.vs+XYBmove XB.vs]', [repmat([1 0 0]', 1, size(YB.vs, 2)) repmat([1 1 0]', 1, size(XB.vs, 2))]', [], dataOutRoot, fileName, 0 );
fileName = 'XY_seg_ori';
fcn_saveUniformSizeModel( [YB.vs XB.vs]', [repmat([1 0 0]', 1, size(YB.vs, 2)) repmat([1 1 0]', 1, size(XB.vs, 2))]', [], dataOutRoot, fileName, 0 );
fileName = 'XY_color';
fcn_saveUniformSizeModel( [YB.vs+XYBmove XB.vs]', [YB.cs XB.cs]', [], dataOutRoot, fileName, 0 );
fileName = 'XY_color_ori';
fcn_saveUniformSizeModel( [YB.vs XB.vs]', [YB.cs XB.cs]', [], dataOutRoot, fileName, 0 );

%%
% Init full set of options %%%%%%%%%%
opt.method='nonrigid'; % use nonrigid registration with lowrank kernel approximation
opt.numeig=30;                 % leave only 30 larges (out of 8171) eigenvectors/values to approximate G
opt.eigfgt=0;                  % use FGT to find the largest eigenvectore/values 

opt.beta=0.05;            % the width of Gaussian kernel (smoothness)
opt.lambda=5;          % regularization weight

opt.viz=1;              % show every iteration
opt.outliers=0.1;       % use 0.7 noise weight
opt.fgt=0;              % use FGT to compute matrix-vector products (2 means to switch to truncated version at the end, see cpd_register)
opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
opt.corresp=0;          % compute correspondence vector at the end of registration (not being estimated by default)

opt.max_it=50;         % max number of iterations
opt.tol=1e-7;           % tolerance
opt.sigma2 = 50;


data.X = X;
data.Y = Y;
data.Xc = Xc;
data.Yc = Yc;
data.Xn = Xn;
data.Yn = Yn;
data.cLambda = 1;
data.nLambda = 1;
opt.auxData = data;
close all;
tic;
[TransformX, ~]=cpd_register(X,Y, opt);
toc;
fmY = TransformX.Y;
diffY = (fmY - Y)';
YBmove = diffY(:, iValidIdx1(Yi));
pclviewer([YB.vs+YBmove XB.vs; repmat([1 0 0]', 1, size(YB.vs, 2)) repmat([1 1 0]', 1, size(XB.vs, 2))]);
pclviewer([YB.vs+YBmove YB.vs; repmat([1 0 0]', 1, size(YB.vs, 2)) repmat([0 1 1]', 1, size(YB.vs, 2))]);
%%


% pclviewer([YB.vs+YBmove]);
pclviewer([YB.vs+YBmove XB.vs; YB.cs XB.cs]);
pclviewer([YB.vs XB.vs; repmat([1 0 0]', 1, size(YB.vs, 2)) repmat([1 1 0]', 1, size(XB.vs, 2))]);

pclviewer([YB.vs+YBmove XB.vs; repmat([1 0 0]', 1, size(YB.vs, 2)) repmat([1 1 0]', 1, size(XB.vs, 2))]);
pclviewer([YB.vs+YBmove YB.vs; repmat([1 0 0]', 1, size(YB.vs, 2)) repmat([1 1 0]', 1, size(YB.vs, 2))]);
pclviewer([XB.vs+XBmove XB.vs; repmat([1 0 0]', 1, size(XB.vs, 2)) repmat([1 1 0]', 1, size(XB.vs, 2))]);

pclviewer([XB.vs; XB.cs]);


dataOutRoot = '/media/mhlee/14A2ACBFA2ACA6A8/E/lab/work/2015_4d_seg/codes/libs/CPD2/examples';
fileName = 'after_move2';
fcn_saveUniformSizeModel( [YB.vs+YBmove XB.vs]', [repmat([1 0 0]', 1, size(YB.vs, 2)) repmat([1 1 0]', 1, size(XB.vs, 2))]', [], dataOutRoot, fileName, 0 );
fileName = 'before_move';
fcn_saveUniformSizeModel( [YB.vs XB.vs]', [repmat([1 0 0]', 1, size(YB.vs, 2)) repmat([1 1 0]', 1, size(XB.vs, 2))]', [], dataOutRoot, fileName, 0 );

pclviewer([config.scene repmat(colors(1,:), size(config.scene, 1), 1) ; fm repmat(colors(2,:), size(config.model, 1), 1) ; config.model repmat(colors(3,:), size(config.model, 1), 1)]');
pclviewer([config.scene repmat(colors(1,:), size(config.scene, 1), 1) ; fm repmat(colors(2,:), size(config.model, 1), 1) ; config.model repmat(colors(4,:), size(config.model, 1), 1)]');
pclviewer([config.scene repmat([1 0 0], size(config.scene, 1), 1) ; config.model repmat([1 1 0], size(config.model, 1), 1)]');
pclviewer([fm repmat([1 0 0], size(config.model, 1), 1) ; config.model repmat([1 1 0], size(config.model, 1), 1)]');
pclviewer([fm repmat([1 0 0], size(config.model, 1), 1) ; config.scene repmat([1 1 0], size(config.scene, 1), 1)]');


pclviewer([XB.vs YB.vs ;XB.cs YB.cs]);

% end
% 
% function [] = saveNonrigidResult()
%     dataOutRoot = '/media/mhlee/14A2ACBFA2ACA6A8/E/lab/work/2015_4d_seg/codes/libs/CPD2/examples';
%     fileName = 'XY_seg';
%     fcn_saveUniformSizeModel( [YB.vs+XYBmove XB.vs]', [repmat([1 0 0]', 1, size(YB.vs, 2)) repmat([1 1 0]', 1, size(XB.vs, 2))]', [], dataOutRoot, fileName, 0 );
%     fileName = 'XY_seg_ori';
%     fcn_saveUniformSizeModel( [YB.vs XB.vs]', [repmat([1 0 0]', 1, size(YB.vs, 2)) repmat([1 1 0]', 1, size(XB.vs, 2))]', [], dataOutRoot, fileName, 0 );
%     fileName = 'XY_color';
%     fcn_saveUniformSizeModel( [YB.vs+XYBmove XB.vs]', [YB.cs XB.cs]', [], dataOutRoot, fileName, 0 );
%     fileName = 'XY_color_ori';
%     fcn_saveUniformSizeModel( [YB.vs XB.vs]', [YB.cs XB.cs]', [], dataOutRoot, fileName, 0 );
% end

