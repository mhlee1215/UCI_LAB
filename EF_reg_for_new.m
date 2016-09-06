% DEMOJRMPCSYNTHETIC   Example of using jrmpc into synthetic data.
%    This example loads numel(theta) views from ./syntheticData/ and calls
%    jrmpc to do the registration. It creates 4 plots one with the initial
%    position of the point sets, one which shows the registration at
%    every iteration, one with the final alignment achieved after maxNumIter
%    iterations and one with the "cleaned up" point sets. Directory 
%    ./syntheticData/ contains 4 partial views from the stanford bunny, 
%    each view is degraded with disparsity noise and outliers. The angles in
%    theta are ground truth angles (same for all 3 axes) used in the 
%    construction.
%
%    $ 18 / 12 / 2014 3:24 PM $

clc
close all
clear all
% g = gpuDevice(1);
% reset(g);

addpath(genpath('libs'));
% addpath(genpath('libs2'));
% run('libs2/vlfeat-0.9.20-bin/vlfeat-0.9.20/toolbox/vl_setup');


% data3DFilePath = '/home/mhlee/Kinect_Logs/2016-02-03.00.klg.ply';
% [tri, pts, data, comments] = ply_read(data3DFilePath, 'tri');
% 
% v = [data.vertex.x data.vertex.y data.vertex.z];
% 
% 
% figure; fv_alligned.Faces = tri;
% fv_alligned.Vertices = v;
% patch(fv_alligned,'FaceColor',       [0.8 0.8 1.0], ...
%          'EdgeColor',       'none',        ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);
% camlight('headlight');
% material('dull');
% axis equal;
% hold on;

%data_set = [1:6 8:14];%[1 2 3 4];
dataRoot = '/home/mhlee/data_from_odroid';
catName = 'LAB_1';
refDataRoot = sprintf('%s/match', dataRoot);
completeDataRoot = sprintf('%s/complete', dataRoot);
refName = sprintf('%s_ref', catName);
refFeature = sprintf('%s.mat', refName);
load(sprintf('%s/%s', refDataRoot, refFeature), 'feature');
refFeaturePath = sprintf('%s/%s', refDataRoot, refFeature);
refFilePath = sprintf('%s/%s.ply', refDataRoot, refName);


fileName = 'LAB_1-2016-07-18_13_14.klg_cvt.ply';
fileName2 = 'LAB_1-2016-07-18_13_16.klg_cvt.ply';
% fileList = dir(sprintf('%s/%s*_cvt.ply', completeDataRoot, catName));

% matchInfo = customXml2struct('/home/mhlee/data_from_odroid/match/match_22_1.mlp');
data_set = {};
data_set.name = {};
data_set.name{1} = sprintf('%s', fileName);
data_set.name{2} = sprintf('%s', fileName2);

d_size = length(data_set.name);
% 
% 
[dataSet, poseSet] = loadEFDataset( completeDataRoot, data_set.name);
dataSetSampled = loadUniformSampling(dataSet, 20);

for i=1:length(dataSetSampled)
    data.vs{i} = dataSetSampled{i}.v;
    data.ns{i} = dataSetSampled{i}.n;
    data.cs{i} = dataSetSampled{i}.c;
    data.is{i} = dataSetSampled{i}.sIndex;
    data.vb{i} = dataSet{i}.v;
    data.nb{i} = dataSet{i}.n;
    data.cb{i} = dataSet{i}.c;
    data.ib{i} = [];
end

data.vs = data.vs';
data.ns = data.ns';
data.cs = data.cs';
data.vb = data.vb';
data.nb = data.nb';
data.cb = data.cb';

featureNew = {};
params_desc.normalRadius=0.2;
params_desc.searchRadius=1.2;
params_desc.searchK = 0;

v1 = data.vs{1}';
c1 = data.cs{1}';
n1 = data.ns{1}';

[k1,~,~] = uniformSubSample(v1, 7, c1);
descriptors = FPFHSDescriptor(v1, c1, k1, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);

featureNew.keypoints = k1;
featureNew.descriptors = descriptors;
featureNew.v = v1;
featureNew.c = c1;
featureNew.n = n1;




% 
% fromIdx = 3;
% toIdx = 1;

v1 = data.vs{1}';%feature.v;
c1 = data.cs{1}';%feature.c;
n1 = data.ns{1}';%feature.n;

v2 = data.vs{2}';%featureNew.v;
c2 = data.cs{2}';%featureNew.n;
n2 = data.ns{2}';%featureNew.c;

params_desc.normalRadius=0.2;
params_desc.searchRadius=1.2;
params_desc.searchK = 0;

[k1,~,~] = uniformSubSample(v1, 5, c1);
f1 = FPFHSDescriptor(v1, c1, k1, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);
% mean(sum(descriptors==0,1))
% k1 = feature.keypoints;
% f1 = feature.descriptors;
 
[k2,~,~] = uniformSubSample(v2, 5, c2);
f2 = FPFHSDescriptor(v2, c2, k2, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);
% k2 = featureNew.keypoints;
% f2 = featureNew.descriptors;

[ initBestR, initBestT ] = featureBasedRegistration( v1, c1, n1, v2, c2, n2, k1, f1, k2, f2, 300 );

% [initBestR, initBestT] = icp(v1', v2', 10000);

v2t = bsxfun(@plus, initBestR*v2', initBestT)';
n2t = bsxfun(@plus, initBestR*n2', initBestT)';
pclviewer([v1 c1 ; v2t c2]');
% pclviewer([v1 c1]');
% pclviewer([v2t c2]');
pclviewer([v1 c1 ; v2 c2]');
% 
return;




V = {};
N = {};
C = {};

V{1} = v1';
C{1} = c1';
N{1} = n1';
V{2} = v2t';
C{2} = c2';
N{2} = n2t';


V = V';
C = C';
N = N';

% V = data.vs(idx);
% N = data.ns(idx);
% C = data.cs(idx);

[R,t,X,S,a,pk,T,TAssigned, TXQ] = joint_align(V,N,C);

params.type = 2;
params.Assigned = TAssigned;
params.K = length(TXQ{1, 1, end});
params.interval = 5;
% params.marker = marker;
% params.markerSize = markerSize;
% params.clrmap = clrmap;
% params.strIdx = strIdx;
params.view = [];%[40 54];
params.pause = 1;
params.TXQ = TXQ;
h = figure;
[TV] = drawTransformation(V, T, params);



vaa = [];
for i = 1:length(TV)
    vaa = [vaa [ TV{i} ;C{i}]];
end
% pclviewer([V{1} ;C{1}]);
pclviewer(vaa);











return;


% idx = 1:d_size;
% 
% V = data.vs(idx);
% N = data.ns(idx);
% C = data.cs(idx);
% 
% M = length(V);
% idx = transpose(1:M);
% strIdx = arrayfun(@(j) sprintf('view %d',j),idx,'uniformoutput',false);
% clrmap = {[1 1 0]; [1 0 1] ; [0 1 1] ; [1 0 0] ; [0 1 0] ; [0 0 1] ; [1 1 1] ; [0 0 0] ; [0.3 0.3 0.3] ; [.5 0.9 0.9] ; [0.8 1 1]; [0 0 0] ...
%     ; [1 1 0]; [1 0 1] ; [0 1 1] ; [1 0 0] ; [0 1 0] ; [0 0 1] ; [1 1 1] ; [0 0 0] ; [0.3 0.3 0.3] ; [.5 0.9 0.9] ; [0.8 1 1]; [0 0 0]};
% clrmap = clrmap(1:M);
% 
% markerSize = {};
% for ii=1:M
%     markerSize{end+1} = 1+ii;
% end
% markerSize = markerSize';
% markerSize = markerSize(1:M);
% 
% marker_set = {'s', 'x', '.', '^'};
% marker = {};
% for i=1:M
%     marker{end+1} = marker_set{mod(i, 4)+1};
% end
% marker = marker';
% 
% % K = 300;
% % randIdx = randperm(size(V{1}', 1));
% % Xin = V{1}(:, randIdx(1:K));
% 
% va = [];
% for i = idx
%     va = [va V{i}];
% end
% [Xin,~,~,~,density]= uniformSubSample( va', 3);
% Xin = Xin(:, find(density > 10));
% K = size(Xin, 2);
% 
% maxNumIter = 100;
% [R,t,X,S,a,pk,T,TAssigned, TXQ] = jrmpc_soft_with_normal(V,Xin, ...
%     'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0);
% 
% % [R,t,X,S,a,pk,T,TAssigned, TXQ] = jrmpc_soft_with_normal(V,Xin, ...
% %     'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0, ...
% %     'normal', N, 'normalLambda', .1, 'color', C, 'colorLambda', .1);
% 
% 
% params.type = 2;
% params.Assigned = TAssigned;
% params.K = K;
% params.interval = 3;
% params.marker = marker;
% params.markerSize = markerSize;
% params.clrmap = clrmap;
% params.strIdx = strIdx;
% params.view = [];%[40 54];
% params.pause = 1;
% params.TXQ = TXQ;
% h = figure;
% [TV] = drawTransformation(V, T, params);
% 
% vaa = [];
% for i = [1 3]%1:length(TV)
%     vaa = [vaa [ TV{i} ;C{i}]];
% end
% % pclviewer([V{1} ;C{1}]);
% pclviewer(vaa);




fromIdx = 3;
toIdx = 1;

v1 = data.vs{toIdx}';
n1 = data.ns{toIdx}';
c1 = data.cs{toIdx}';

v2 = data.vs{fromIdx}';
n2 = data.ns{fromIdx}';
c2 = data.cs{fromIdx}';

params_desc.normalRadius=0.2;
params_desc.searchRadius=1.2;
params_desc.searchK = 0;

[k1,~,~] = uniformSubSample(v1, 5, c1);
descriptors = FPFHSDescriptor(v1, c1, k1, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);
mean(sum(descriptors==0,1))
 
[k2,~,~] = uniformSubSample(v2, 5, c2);
f2 = FPFHSDescriptor(v2, c2, k2, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);


[ initBestR, initBestT ] = featureBasedRegistration( v1, c1, n1, v2, c2, n2, k1, descriptors, k2, f2, 300 );
v_m2_best11 = bsxfun(@plus, initBestR*v2', initBestT)';
pclviewer([v1 c1 ; v_m2_best11 c2]');



