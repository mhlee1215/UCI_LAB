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
dataRoot = '/home/mhlee/data_from_odroid/complete';

catName = 'LAB_1';
fileList = dir(sprintf('%s/%s*_cvt.ply', dataRoot, catName));

% matchInfo = customXml2struct('/home/mhlee/data_from_odroid/match/match_22_1.mlp');
data_set = {};
data_set.name = {};
for i=1:3:length(fileList);
    data_set.name{end+1} = fileList(i).name;
end
% % data_set.R = {};
% % data_set.t = {};
% for i=1:length(matchInfo.MeshLabProject.MeshGroup.MLMesh)
%     data_set.name{i} = matchInfo.MeshLabProject.MeshGroup.MLMesh{i}.Attributes.filename;
%     MatText = matchInfo.MeshLabProject.MeshGroup.MLMesh{i}.MLMatrix44.Text;
%     
%     curR = [];
%     MatRows = strsplit(strtrim(MatText), '\n');
%     for j = 1:length(MatRows)
%         elements = strsplit(strtrim(MatRows{j}), ' ');
%         for k = 1:length(elements)
%             curR(end+1) = str2num(elements{k});
%         end
%     end
%     curR = reshape(curR, 4, 4)';
%     data_set.R{i} = curR(1:3, 1:3);
%     data_set.t{i} = curR(1:3, 4);
% end

% data_set = {'LAB_1-2016-07-22_12_32.klg_cvt', 'LAB_2-2016-07-22_12_41.klg_cvt'};

d_size = length(data_set.name);
% 
% 
[dataSet, poseSet] = loadEFDataset( dataRoot, data_set.name);
dataSetSampled = loadUniformSampling(dataSet, 10);
dataSetSampled2 = loadUniformSampling(dataSet, 20);
% 


% V = {};
% N = {};
% C = {};
% Vall = {};
% Nall = {};
% Call = {};
for i=1:length(dataSetSampled)
    data.vs{i} = dataSetSampled{i}.v;%bsxfun(@plus, data_set.R{i}*dataSetSampled{i}.v, data_set.t{i});
    data.ns{i} = dataSetSampled{i}.n;
    data.cs{i} = dataSetSampled{i}.c;
    data.is{i} = dataSetSampled{i}.sIndex;
%     data.vm{i} = dataSetSampled2{i}.v;
%     data.nm{i} = dataSetSampled2{i}.n;
%     data.cm{i} = dataSetSampled2{i}.c;
%     data.im{i} = dataSetSampled2{i}.sIndex;
    data.vb{i} = dataSet{i}.v;
    data.nb{i} = dataSet{i}.n;
    data.cb{i} = dataSet{i}.c;
    data.ib{i} = [];
end

data.vs = data.vs';
data.ns = data.ns';
data.cs = data.cs';
% data.vm = data.vm';
% data.nm = data.nm';
% data.cm = data.cm';
data.vb = data.vb';
data.nb = data.nb';
data.cb = data.cb';

feature = {};
params_desc.normalRadius=0.2;
params_desc.searchRadius=1.2;
params_desc.searchK = 0;

for i=1:d_size
    v1 = data.vs{i}';
    c1 = data.cs{i}';
    n1 = data.ns{i}';
    
    [keypoints,~,~] = uniformSubSample(v1, 5, c1);
    descriptors = FPFHSDescriptor(v1, c1, keypoints, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);
    
    feature{i}.keypoints = keypoints;
    feature{i}.descriptors = descriptors;
end

% return;


idx = 1:d_size;

V = data.vs(idx);
N = data.ns(idx);
C = data.cs(idx);

M = length(V);
idx = transpose(1:M);
strIdx = arrayfun(@(j) sprintf('view %d',j),idx,'uniformoutput',false);
clrmap = {[1 1 0]; [1 0 1] ; [0 1 1] ; [1 0 0] ; [0 1 0] ; [0 0 1] ; [1 1 1] ; [0 0 0] ; [0.3 0.3 0.3] ; [.5 0.9 0.9] ; [0.8 1 1]; [0 0 0] ...
    ; [1 1 0]; [1 0 1] ; [0 1 1] ; [1 0 0] ; [0 1 0] ; [0 0 1] ; [1 1 1] ; [0 0 0] ; [0.3 0.3 0.3] ; [.5 0.9 0.9] ; [0.8 1 1]; [0 0 0]};
clrmap = clrmap(1:M);

markerSize = {};
for ii=1:M
    markerSize{end+1} = 1+ii;
end
markerSize = markerSize';
markerSize = markerSize(1:M);

marker_set = {'s', 'x', '.', '^'};
marker = {};
for i=1:M
    marker{end+1} = marker_set{mod(i, 4)+1};
end
marker = marker';

% K = 300;
% randIdx = randperm(size(V{1}', 1));
% Xin = V{1}(:, randIdx(1:K));

va = [];
for i = idx
    va = [va V{i}];
end
[Xin,~,~,~,density]= uniformSubSample( va', 3);
Xin = Xin(:, find(density > 10));
K = size(Xin, 2);

maxNumIter = 100;
[R,t,X,S,a,pk,T,TAssigned, TXQ] = jrmpc_soft_with_normal(V,Xin, ...
    'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0);

% [R,t,X,S,a,pk,T,TAssigned, TXQ] = jrmpc_soft_with_normal(V,Xin, ...
%     'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0, ...
%     'normal', N, 'normalLambda', .1, 'color', C, 'colorLambda', .1);


params.type = 2;
params.Assigned = TAssigned;
params.K = K;
params.interval = 3;
params.marker = marker;
params.markerSize = markerSize;
params.clrmap = clrmap;
params.strIdx = strIdx;
params.view = [];%[40 54];
params.pause = 1;
params.TXQ = TXQ;
h = figure;
[TV] = drawTransformation(V, T, params);

vaa = [];
for i = [1 3]%1:length(TV)
    vaa = [vaa [ TV{i} ;C{i}]];
end
% pclviewer([V{1} ;C{1}]);
pclviewer(vaa);




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

[keypoints,~,~] = uniformSubSample(v1, 5, c1);
descriptors = FPFHSDescriptor(v1, c1, keypoints, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);
mean(sum(descriptors==0,1))
 
[k2,~,~] = uniformSubSample(v2, 5, c2);
f2 = FPFHSDescriptor(v2, c2, k2, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);


[ initBestR, initBestT ] = featureBasedRegistration( v1, c1, n1, v2, c2, n2, keypoints, descriptors, k2, f2, 300 );
v_m2_best11 = bsxfun(@plus, initBestR*v2', initBestT)';
pclviewer([v1 c1 ; v_m2_best11 c2]');



