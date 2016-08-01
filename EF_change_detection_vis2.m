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
addpath(genpath('libs2'));
run('libs2/vlfeat-0.9.20-bin/vlfeat-0.9.20/toolbox/vl_setup');


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
dataRoot = '/home/mhlee/Kinect_Logs_test';
% data_set = {'2016-02-09.01','2016-02-09.02','2016-02-09.03'};
% data_set = {'2016-04-29.03', '2016-04-29.05', '2016-04-29.06', '2016-04-29.00','2016-04-29.08'};
%Wrong init pos : '2016-04-29.06'
% data_set = {'2016-04-29.00', '2016-04-29.01', '2016-04-29.02', '2016-04-29.03', '2016-04-29.05', '2016-04-29.08'};%, '2016-04-29.05'};%, '2016-04-29.00','2016-04-29.08'};
data_set = {'2016-04-29.00', '2016-04-29.01', '2016-04-29.02', '2016-04-29.03', '2016-04-29.05', '2016-04-29.08'};%, '2016-04-29.05'};%, '2016-04-29.00','2016-04-29.08'};
d_size = length(data_set);


[dataSet, poseSet] = loadEFDataset( dataRoot, data_set);
dataSetSampled = loadUniformSampling(dataSet, 20);

V = {};
N = {};
C = {};
for i=1:length(dataSetSampled)
    V{i} = dataSetSampled{i}.v;
    N{i} = dataSetSampled{i}.n;
    C{i} = dataSetSampled{i}.c;
end
V = V';
N = N';
C = C';


resultPath = 'results/EF_vis3';
mkdir(resultPath);


maxNumIter = 100;                    
                     
% latent angles, rotating V{j} by theta(j) reprojects it to V{1} rotated by
% theta(1). used to quantify the accuracy of the estimated R
% theta = [0; pi/20; pi/10; pi/6];
% theta = theta(1:d_size);

% number of views, M files must be found in the directory ./syntheticData 
M = d_size;%numel(theta);

% cell with indexes 1:M, used as suffixes on view's filenames
idx = transpose(1:M);

% string-labels for legends in subsequent plots
strIdx = arrayfun(@(j) sprintf('view %d',j),idx,'uniformoutput',false);

% fprintf(' Data loading...\n');

% load the views, file view<j>.txt corresponds to theta(j)
% cutView*.txt is a partial view as described in the paper, while view*.txt
% "sees" the whole surface (again downsampled and noisy)

%V = arrayfun(@(j) dlmread(sprintf('./syntheticData/view%d.txt',j),' ')',idx,'uniformoutput',false);
% V = arrayfun(@(j) dlmread(sprintf('libs/JRMPC_v0.9.4/syntheticData/cutView%d.txt',j),' ')',idx,'uniformoutput',false);

% sampleNum = 5000;



%-0.55 to 0.55, -0.44 to 0.44


% return;


% ground truth rotation matrices Rgt{j}*V{j} is aligned with Rgt{1}*R{1}
% Rgt = arrayfun(@(theta) angle2rotation(theta),theta,'uniformoutput',false);

% colors for each view
% clrmap = {[1 .1412 0]; [.1373 .4196 .5569]; [0 0 1]; [.8039 .6078 .1137]};
clrmap = {[1 1 0]; [1 0 1] ; [0 1 1] ; [1 0 0] ; [0 1 0] ; [0 0 1] ; [1 1 1] ; [0 0 0] ; [0.3 0.3 0.3] ; [.5 0.9 0.9] ; [0.8 1 1]; [0 0 0] ...
    ; [1 1 0]; [1 0 1] ; [0 1 1] ; [1 0 0] ; [0 1 0] ; [0 0 1] ; [1 1 1] ; [0 0 0] ; [0.3 0.3 0.3] ; [.5 0.9 0.9] ; [0.8 1 1]; [0 0 0]};
clrmap = clrmap(1:d_size);

% markerSizes 
% markerSize = {7; 70; 12; 10};
% markerSize = {2; 3; 4; 5};
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
    marker{end+1} = marker_set{mod(i, 4)+1};
end
marker = marker';
% marker = {'s'; 'x'; '.'; '^';'s'; 'x'; '.'; '^';'s'; 'x'; '.'; '^';'s'...
%     ; 'x'; '.'; '^';'s'; 'x'; '.'; '^';'s'; 'x'; '.'; '^'};
% marker = marker(1:d_size);

% initialize GMM means Xin, using random sampling of a unit sphere. Choose
% your own initialization. You may want to initialize Xin with some of the
% sets.

% set K as the 50% of the median cardinality of the views
K = ceil(0.5*median(cellfun(@(V) size(V,2),V))); 
K = 1000;

rGMM_color = bsxfun(@plus, rand(K, 3)./2, [0.3 0.3 0.3]);

% sample the unit sphere, by randomly selecting azimuth / elevation angles
az = 2*pi*rand(1,K);
el = 2*pi*rand(1,K);

%points on a unit sphere
% Xin = [cos(az).*cos(el); sin(el); sin(az).*cos(el)];% (unit) polar to cartesian conversion
% 
% Xin = Xin/10; % it is good for the initialization to have initial cluster centers at the same order with the points
% % since sigma is automatically initialized based on X and V
% 
% meanData = cellfun(@(a) mean(a')',V,'uniformoutput',false); % 1 x K rows
% meanData = mean(cat(3,meanData{:}),3);
% % Xin = Xin ./ repmat(sqrt(var(Xin'))', 1, size(Xin, 2));
% Xin = Xin + repmat(meanData, 1, size(Xin, 2));

randIdx = randperm(size(V{1}', 1));
% K = 1000;
Xin = V{1}(:, randIdx(1:K));


% show the initial position of the point clouds
h=figure;
hold on,grid on

% make the legend
title('Initial position of the point clouds','fontweight','bold','fontsize',12);

hg1 = cellfun(@(V,clrmap,marker,markerSize) scatter3(V(1,:),V(2,:),V(3,:),markerSize.*1.5,[0 0 0],'filled'),V,clrmap,marker,markerSize, 'UniformOutput', false);
hg1 = cellfun(@(V,clrmap,marker,markerSize) scatter3(V(1,:),V(2,:),V(3,:),markerSize,clrmap,'filled'),V,clrmap,marker,markerSize, 'UniformOutput', false);

legend(strIdx{:});

set(1,'position',get(1,'position')+[-260 0 0 0]);

% set(gca,'fontweight','bold','children',hg1);
set(gca,'fontweight','bold');

view([40 54])
scatter3(Xin(1,:),Xin(2,:),Xin(3,:),'k')
hold off; drawnow

saveas(h, sprintf('%s/init_pos', resultPath));
close(h);

fprintf('Data registration... \n\n');

tic;
% call JRMPC (type jrmpc with no arguments to see the documentation).
% [R,t,X,S,a,pk,T] = jrmpc(V,Xin,'maxNumIter',maxNumIter,'gamma',0.1);
% [R,t,X,S,a,pk,T] = jrmpc(V,Xin,'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 1);

% tic;
% [R,t,X,S,a,pk,T] = jrmpc2_gpu(V,Xin,'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0);
% toc;
% tic
% [R,t,X,S,a,pk,T,TAssigned, TXQ] = jrmpc_soft(V,Xin,'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0);
% toc


maxNumIter = 100;
% [R,t,X,S,a,pk,T,TAssigned, TXQ] = jrmpc_soft_with_normal(V,Xin, ...
%     'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0, ...
%     'normal', N, 'normalLambda', .2, 'color', C, 'colorLambda', .0);

%Parameter Test
% jrResult = {};
% for i=16:18
%     normalLambda = 0.9 - 0.05*i;
%     jrResult{i}.normalLambda = normalLambda;
%     [R,t,X,S,a,pk,T,TAssigned, TXQ] = jrmpc_soft_with_normal(V,Xin, ...
%         'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0, ...
%         'normal', N, 'normalLambda', normalLambda);
%     jrResult{i}.TAssigned = TAssigned;
% end
% 
GMM_color = bsxfun(@plus, rand(params.K, 3)./2, [0.3 0.3 0.3]);
pcl_model = [];
ii=10;
for i=1:M
    clustAssin = cell2mat(jrResult{ii}.TAssigned(i,end));
    pcl_model = [pcl_model [TV{i} ;GMM_color(clustAssin, :)']];
end
pclviewer(pcl_model);



[R,t,X,S,a,pk,T,TAssigned, TXQ] = jrmpc_soft_with_normal(V,Xin, ...
    'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0, ...
    'normal', N, 'normalLambda', 0.15);


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
[TV] = drawTransformation(V, T, params);


% GMM_color = bsxfun(@plus, rand(params.K, 3)./2, [0.3 0.3 0.3]);
% pcl_model = [];
% for i=1:M
%     clustAssin = cell2mat(TAssigned(i,end));
%     pcl_model = [pcl_model [TV{i} ;GMM_color(clustAssin, :)']];
% end
% pclviewer(pcl_model);
% 
% 
% [R0,t0,X0,S0,a0,pk0,T0,TAssigned0, TXQ0] = jrmpc_soft_with_normal(V,Xin, ...
%     'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0);
% pcl_model = [];
% for i=1:M
%     clustAssin = cell2mat(TAssigned0(i,end));
%     pcl_model = [pcl_model [TV{i} ;GMM_color(clustAssin, :)']];
% end
% pclviewer(pcl_model);

return;

% toc;
% tic;
% [R,t,X,S,a,pk,T] = jrmpc2_gpu(V,Xin,'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 1);
% toc;

% measure and display convergency, view 1 is ommited as is the referential.
% fprintf('                  ||Rgt{j} - R{j}^T*R{1}||_F                  \n');
% 
% fprintf('______________________________________________________________\n');
% 
% fprintf('Set  :'),for j=2:M,fprintf('    %d    ',j),end,fprintf('\n');
% 
% % fprintf('Error:'),for j=2:M,fprintf('  %.4f ',norm(Rgt{j}-R{j}'*R{1},'fro'));end
% 
% fprintf('\n');


% visualize the registration process, see documentation of jrmpc for T.

iter=maxNumIter;
TV_all =       cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V_all',T(:,1,iter),T(:,2,iter),'uniformoutput',false);

TV_small_all = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V,T(:,1,iter),T(:,2,iter),'uniformoutput',false);


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
[TV] = drawTransformation(V, T, params);
axis equal;



return;

GMM_color = bsxfun(@plus, rand(params.K, 3)./2, [0.3 0.3 0.3]);
pcl_model = [];
for i=1:M
    clustAssin = cell2mat(TAssigned(i,end));
    pcl_model = [pcl_model [TV{i} ;GMM_color(clustAssin, :)']];
end
pclviewer(pcl_model);




% srcId = 1;
% targetId = 1;
% 
% vertices = V{targetId};
% pSet = Pose_Set{srcId};
% 
% visible = getCameraVisiblePoints(vertices, pSet);
% 
% colors = C2{targetId};
% scatter3(vertices(1, find(visible)), vertices(2, find(visible)), vertices(3, find(visible)), 8, colors(:, find(visible))', 'filled');
% view(182, -53)
% axis([-5 5 -3 3 -5 5]);
% 
% scatter3(vertices(1, find(visible==0)), vertices(2, find(visible==0)), vertices(3, find(visible==0)), 8, colors(:, find(visible==0))', 'filled');
% view(182, -53)
% axis([-5 5 -3 3 -5 5]);

% refineIterMax = 30;
% [R2,t2,X2,S2,a2,pk2,T2,TAssigned2, TXQ2] = jrmpc_soft(V,X,'maxNumIter',refineIterMax, 'updatepriors', 0, 'r', R, 't', t, 's', S, 'initialpriors', pk', 'updateVis', 0);

% observedAll = {};
unObservedAllX = zeros(M, K);
observedAllX = zeros(M, K);
observedAllX2 = zeros(M, K);
observedAllXBack = zeros(M, K);
existAllX = zeros(M, K);

for srcId = 1:M
    srcId
% vertices = [TV{1} TV{2} TV{3} TV{4}];
vertices = V{srcId};%[TV{1} TV{2} TV{3} TV{4}];
pSet = poseSet{srcId};
paramsVis.R = cell2mat(T(srcId,1,end));
paramsVis.T = cell2mat(T(srcId,2,end));
colors = dataSetSampled{srcId}.c;%[C2{1} C2{2} C2{3} C2{4}];
clustAssign = TAssigned(:,end);

paramsVis.srcId = srcId;
paramsVis.colors = colors;
paramsVis.visParam = 1;
paramsVis.camInterval = 30;
paramsVis.X = X;
paramsVis.TV = TV;
paramsVis.clustAssign = clustAssign;

[observed, observedX2, existX, visibleXBack] = getCameraVisiblePoints2(vertices, pSet, paramsVis);


observedX = zeros(size(X, 2), 1);
for i=1:size(X, 2)
    pIdx = find(clustAssign{srcId} == i);
    if sum(observed(pIdx)) / length(pIdx) > 0.7
        observedX(i) = 1;
    end
end

% pclviewer(X(:,find(observedX)));
% 
% scatter3(vertices(1, find(observed)), vertices(2, find(observed)), vertices(3, find(observed)), 8, colors(:, find(observed))', 'filled');
% view(182, -53)
% axis([-5 5 -3 3 -5 5]);
% 
% scatter3(vertices(1, find(observed==0)), vertices(2, find(observed==0)), vertices(3, find(observed==0)), 8, colors(:, find(observed==0))', 'filled');
% view(182, -53)
% axis([-5 5 -3 pcl3 -5 5]);

% observedAll{srcId} = observed';
observedAllX(srcId,:)= observedX';
observedAllX2(srcId,:) = observedX2';
observedAllXBack(srcId, :) = visibleXBack';
existAllX(srcId,:) = existX';
unObservedAllX(srcId,:) = min(observedX+existX, 1)';

end



% return;

% pclviewer([TV{3} ;SEG_color(observedAll{3}+1,:)']);


SEG_color = [[0.3 0.3 0.3];[0.9 0.1 0.1]; [0.1 0.9 0.1] ; [0.1 0.1 0.9]];


%view unobserved
src = 1;
pcl_model = [];
for i=1:M
    clustAssin = cell2mat(TAssigned(i,end));
    unObservedVec = ~observedAllX(src, clustAssin);%max(observedAllX(src,clustAssin), observedAllX2(src,clustAssin));
    existVec = existAllX(src, clustAssin);%max(observedAllX(src,clustAssin), observedAllX2(src,clustAssin));
%     observedVec = max(~unObservedVec - unExistVec, 0);
    visibleX = observedAllX(src, clustAssin);
    existX = existAllX(src, clustAssin);
    visibleXBack = observedAllXBack(src, clustAssin);
    catVec = ones(size(existVec));
    catVec(find(~existX)) = 2;
    catVec(find(visibleXBack)) = 3;
%     catVec(find(observedVec)) = 4;
   
    pcl_model = [pcl_model [TV{i} ;SEG_color(catVec, :)']];
end
pclviewer(pcl_model);

for src = 1%:M
    h=figure; hold on;
    for i=1:M       
        scatter3(TV{i}(1,:)', TV{i}(2,:)', TV{i}(3,:)', 8, C2{i}', 'filled');
    end
    title(sprintf('colorAll viewpoint:%d', src));
    axis equal;
    view(-184, -27);
    axis([-4.5 3 -2.5 2 -4 3]);
    saveas(h, sprintf('%s/colorAll_%d_1.png', resultPath, src));
    
    view(177, -52);
    axis([-4.5 3 -2.5 2 -4 3]);
    saveas(h, sprintf('%s/colorAll_%d_2.png', resultPath, src));
    close(h);
end

for src = 1:M
    h=figure; hold on;
    i = src;
%     for i=1:M       
        scatter3(TV{i}(1,:)', TV{i}(2,:)', TV{i}(3,:)', 8, C2{i}', 'filled');
%     end
    title(sprintf('color viewpoint:%d', src));
    axis equal;
    view(-184, -27);
    axis([-4.5 3 -2.5 2 -4 3]);
    saveas(h, sprintf('%s/color_%d_1.png', resultPath, src));
    
    view(177, -52);
    axis([-4.5 3 -2.5 2 -4 3]);
    saveas(h, sprintf('%s/color_%d_2.png', resultPath, src));
    close(h);
end

for src = 1:M
    h=figure; hold on;
    for i=1:M
        clustAssin = cell2mat(TAssigned(i,end));
        observedVec = unObservedAllX(src, clustAssin);%max(observedAllX(src,clustAssin), observedAllX2(src,clustAssin));
        scatter3(TV{i}(1,:)', TV{i}(2,:)', TV{i}(3,:)', 8, SEG_color(observedVec+1, :), 'filled');
    end
    title(sprintf('Unobserved from viewpoint:%d', src));
    axis equal;
    view(-184, -27);
    axis([-4.5 3 -2.5 2 -4 3]);
    saveas(h, sprintf('%s/unobserved_%d_1.png', resultPath, src));
    
    view(177, -52);
    axis([-4.5 3 -2.5 2 -4 3]);
    saveas(h, sprintf('%s/unobserved_%d_2.png', resultPath, src));
    close(h);
end

for src = 1:M
    h=figure; hold on;
    for i=1:M
        clustAssin = cell2mat(TAssigned(i,end));
        observedVec = existAllX(src, clustAssin);%max(observedAllX(src,clustAssin), observedAllX2(src,clustAssin));
        scatter3(TV{i}(1,:)', TV{i}(2,:)', TV{i}(3,:)', 8, SEG_color(observedVec+1, :), 'filled');
    end
    title(sprintf('Exist from viewpoint:%d', src));
    axis equal;
    view(-184, -27);
    axis([-4.5 3 -2.5 2 -4 3]);
    saveas(h, sprintf('%s/exist_%d_1.png', resultPath, src));
    
    view(177, -52);
    axis([-4.5 3 -2.5 2 -4 3]);
    saveas(h, sprintf('%s/exist_%d_2.png', resultPath, src));
    close(h);
end

%view exist
src = 1;
pcl_model = [];
for i=1:M
    clustAssin = cell2mat(TAssigned(i,end));
    observedVec = existAllX(src,clustAssin);
    pcl_model = [pcl_model [TV{i} ;SEG_color(observedVec+1, :)']];
end
pclviewer(pcl_model);

%view observed
src = 3;
pcl_model = [];
for i=1:M
    clustAssin = cell2mat(TAssigned(i,end));
    observedVec = observedAllX(src,clustAssin);
    pcl_model = [pcl_model [TV{i} ;SEG_color(observedVec+1, :)']];
end
pclviewer(pcl_model);

%color view
pcl_model = [];
for i=4
    pcl_model = [pcl_model [TV{i} ; dataSetSampled{i}.c]];
end
pclviewer(pcl_model);



segId = 1;

[ segmentedPoints, coloredCloud ] = Segmentation(V{segId}', C2{segId}', 15, '');
colorSet = unique(coloredCloud(4:6,:)', 'rows')';
% pclviewer(coloredCloud);
s=cellfun(@size,segmentedPoints,'uniform',false);
[trash is]=sortrows(cat(1,s{:}),-[1 2]);
segmentSorted = segmentedPoints(is);
segmentIndex = getIndexFromVertices(V{segId}', segmentSorted);

%Get major clusters
clustAssin = cell2mat(TAssigned(segId,end));
clusterList = clustAssin(segmentIndex{1});

clusterSum = accumarray(clusterList', ones(1,length(clusterList))', [], @sum);
clusterSumAll = zeros(max(clusterList'),1);
clusterSumAll(1:length(clusterSum)) = clusterSum;

sortedCluster = sortrows([(1:max(clusterList))' clusterSumAll], 2);
sortedCluster(end:-1:1) = sortedCluster;
sortedClusterCrop = sortedCluster(1:min(find(sortedCluster(:,1)==0))-1,:);
majorCluster = sortedClusterCrop(1,1);









fv2 = fv;
for i=1:d_size
    fv2{i}.Vertices = TV_all{i}';
end


refineIterMax = 30;
[R2,t2,X2,S2,a2,pk2,T2,TAssigned2, TXQ2] = jrmpc_soft(V,X,'maxNumIter',refineIterMax, 'updatepriors', 0, 'r', R, 't', t, 's', S, 'initialpriors', pk', 'updateVis', 0);

ajk2 = cellfun(@(a) sum(a), a2, 'uniformoutput', false);
ajkMat2 = cell2mat(ajk2);


src = 4;
[ segmentedPoints, coloredCloud ] = Segmentation(V{src}', C2{src}', 100, '');
[clustCent,point2cluster,clustMembsCell] = ...
        MeanShiftCluster([ajkMat2;X2.*30], 50);
GMM_color = bsxfun(@plus, rand(params.K, 3)./2, [0.3 0.3 0.3]);
clustAssign = cell2mat(TAssigned2(src,end));
pclviewer([TV{src} ; GMM_color(point2cluster(clustAssign), :)']);
pclviewer(coloredCloud);


















return;

%function point = projPointOnLine3d(point, line)
%function d = distancePointLine3d(point, line)








pclviewer([TV{4} ; C2{4}])
pclviewer(X(:,:));
pclviewer(X(:,find(~observedAllX(4,:))));



fv2 = fv;
for i=1:d_size
    fv2{i}.Vertices = TV_all{i}';
end


id = 3;
clear 'paramsPlotGMMK';
plotType = 2;
paramsPlotGMMK.setId = id;
paramsPlotGMMK.gmmk = [];
paramsPlotGMMK.assignments = [];
paramsPlotGMMK.maxIdx = TAssigned;
paramsPlotGMMK.fv = fv2;
paramsPlotGMMK.I2 = I2;
paramsPlotGMMK.type = plotType;
paramsPlotGMMK.gmmk_gmmidx = [];

h=plotGMMKv2( paramsPlotGMMK );
view(179, -36)


% [visiblePtInds, inVisiblePtInds, inVisibleVector]=HPR(p,C,param)

% return;
% ajk = cellfun(@(a) sum(a), a, 'uniformoutput', false);
% ajkMat = cell2mat(ajk);
% 
% epsilonST = [0.0005];
% unObservedAll = ~observedAll;
% 
% ajkMat2 = ajkMat + unObservedAll.*30;
% 
% [ st ] = genUniformDist( ajkMat, epsilonST );
% vis = genVisFromPeriod(st, M, K, epsilonST);
% vis = cell2mat(vis');
% 
% id = 2;
% plotType = 1;
% paramsPlotGMMK.setId = id;
% paramsPlotGMMK.gmmk = [];
% paramsPlotGMMK.assignments = [];
% paramsPlotGMMK.maxIdx = TAssigned;
% paramsPlotGMMK.fv = fv2;
% paramsPlotGMMK.I2 = I2;
% paramsPlotGMMK.type = plotType;
% paramsPlotGMMK.gmmk_gmmidx = find( (st(:,2)-st(:,1)) == M-1 );
% 
% h=plotGMMKv2( paramsPlotGMMK );
% view(179, -36)
% title(sprintf('Without vis Term, setID=%d, e=%.4f', id, epsilonST));


fv2 = fv;
for i=1:d_size
    fv2{i}.Vertices = TV_all{i}';
end


refineIterMax = 30;
[R2,t2,X2,S2,a2,pk2,T2,TAssigned2, TXQ2] = jrmpc_soft(V,X,'maxNumIter',refineIterMax, 'updatepriors', 0, 'r', R, 't', t, 's', S, 'initialpriors', pk', 'updateVis', 0);

ajk2 = cellfun(@(a) sum(a), a2, 'uniformoutput', false);
ajkMat2 = cell2mat(ajk2);
ajkMat3 = ajkMat2 + unObservedAllX.*100;%repmat(mean(ajkMat2), M, 1);

[ st2 ] = genUniformDist( ajkMat3, epsilonST );
% vis2 = genVisFromPeriod(st2, M, K, epsilonST);
% vis2 = cell2mat(vis2');

id = 3;
clear 'paramsPlotGMMK';
plotType = 2;
paramsPlotGMMK.setId = id;
paramsPlotGMMK.gmmk = [];
paramsPlotGMMK.assignments = [];
paramsPlotGMMK.maxIdx = TAssigned2;
paramsPlotGMMK.fv = fv2;
paramsPlotGMMK.I2 = I2;
paramsPlotGMMK.type = plotType;
paramsPlotGMMK.gmmk_gmmidx = find( (st2(:,2)-st2(:,1)) == M-1 );

h=plotGMMKv2( paramsPlotGMMK );
view(179, -36)



return;



refineIterMax = 30;
[R2,t2,X2,S2,a2,pk2,T2,TAssigned2, TXQ2] = jrmpc_soft(V,X,'maxNumIter',refineIterMax, 'updatepriors', 0, 'r', R, 't', t, 's', S, 'initialpriors', pk', 'updateVis', 0);

fv2 = fv;
for i=1:d_size
    fv2{i}.Vertices = TV_all{i}';
end

pclviewer([fv2{1}.Vertices fv2{1}.FaceVertexCData ; fv2{2}.Vertices fv2{2}.FaceVertexCData ; fv2{3}.Vertices fv2{3}.FaceVertexCData ; fv2{4}.Vertices fv2{4}.FaceVertexCData]')

% return;

%Draw changed region
for epsilonST = [0.0005]
    
    
    
    % [R2,t2,X2,S2,a2,pk2,T2] = jrmpc_vis(V,Xin,'maxNumIter',refineIterMax,'gamma',0.1, 'updatepriors', 0, 'r', R, 't', t, 's', S, 'vis', A_binary_c);


    % [R3,t3,X3,S3,a3,pk3,T3,TAssigned3, TXQ3] = jrmpc_soft(V,X,'maxNumIter',refineIterMax, 'updatepriors', 0, 'r', R, 't', t, 's', S, 'initialpriors', pk', 'updateVis', 0, 'vis', A_binary_c);
    % [R4,t4,X4,S4,a4,pk4,T4,TAssigned4, TXQ4] = jrmpc_soft(V,X,'maxNumIter',refineIterMax, 'updatepriors', 0, 'r', R, 't', t, 's', S, 'initialpriors', pk', 'updateTR', 0);
    [R5,t5,X5,S5,a5,pk5,T5,TAssigned5, TXQ5] = jrmpc_soft(V,X,'maxNumIter',refineIterMax, 'updatepriors', 0, 'r', R, 't', t, 's', S, 'initialpriors', pk', 'epsilonST', epsilonST);


    %Without vis term, just 30 more refine with similar setting

    ajk2 = cellfun(@(a) sum(a), a2, 'uniformoutput', false);
    ajkMat2 = cell2mat(ajk2);

    [ st2 ] = genUniformDist( ajkMat2, epsilonST );
    vis2 = genVisFromPeriod(st2, M, K, epsilonST);
    vis2 = cell2mat(vis2');

    for id=1:d_size
        clear 'paramsPlotGMMK';
        plotType = 1;
        paramsPlotGMMK.setId = id;
        paramsPlotGMMK.gmmk = [];
        paramsPlotGMMK.assignments = [];
        paramsPlotGMMK.maxIdx = TAssigned2;
        paramsPlotGMMK.fv = fv2;
        paramsPlotGMMK.I2 = I2;
        paramsPlotGMMK.type = plotType;
        paramsPlotGMMK.gmmk_gmmidx = find( (st2(:,2)-st2(:,1)) == M-1 );
        
        h=plotGMMKv2( paramsPlotGMMK );
        view(179, -36)
        title(sprintf('Without vis Term, setID=%d, e=%.4f', id, epsilonST));
        saveas(h, sprintf('%s/WO_%d_%d_%.4f.png', resultPath, id, plotType, epsilonST));
        close(h);
    end

    %With vis term, just 30 more refine with similar setting
    % epsilonST = 0.05;
    ajk5 = cellfun(@(a) sum(a), a5, 'uniformoutput', false);
    ajkMat5 = cell2mat(ajk5);

    [ st5 ] = genUniformDist( ajkMat5, epsilonST );
    vis5 = genVisFromPeriod(st5, M, K, epsilonST);
    vis5 = cell2mat(vis5');

    for id=1:d_size
        clear 'paramsPlotGMMK';
        plotType = 2;
        paramsPlotGMMK.setId = id;
        paramsPlotGMMK.gmmk = [];
        paramsPlotGMMK.assignments = [];
        paramsPlotGMMK.maxIdx = TAssigned5;
        paramsPlotGMMK.fv = fv2;
        paramsPlotGMMK.I2 = I2;
        paramsPlotGMMK.type = plotType;
        paramsPlotGMMK.gmmk_gmmidx = find( (st5(:,2)-st5(:,1)) == M-1 );
 
        h=plotGMMKv2( paramsPlotGMMK );
        view(179, -36)
        title(sprintf('With vis Term, setID=%d, e=%.4f', id, epsilonST));
        saveas(h, sprintf('%s/W_%d_%d_%.4f.png', resultPath, id, plotType, epsilonST));
        close(h);
    end

end

return;



for DataSetId = 1:14
h=figure;
fv_alligned.Faces = fv{DataSetId}.small.faces;
fv_alligned.Vertices = bsxfun(@plus, T{DataSetId,1,end}*fv{DataSetId}.small.vertices', T{DataSetId,2,end})';
patch(fv_alligned,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');
axis equal;
hold on;


ocl = zeros(M, K);
tol = 0.1;
nHit = 0;
for j=DataSetId%:M    
    RMat = T{j,1,end};
    TMat = T{j,2,end};
    fv_small = fv{j}.small;
    fv_small.vertices = bsxfun(@plus, T{j,1,end}*fv_small.vertices', T{j,2,end})';
    isObserved = ones(1, K).*2;     %Initialize as 2 (unknown)
    
    curPoseSet = Pose_Set{j};
    
    
    poseMat = cell2mat(Pose_Set{j})';
    poseMat = bsxfun(@plus, RMat*poseMat', TMat)';
    scatter3(poseMat(:,1), poseMat(:,2), poseMat(:,3), 50, 'm', 'filled')
    [poseRep] = kmeans(poseMat', 30);
    poseRep = poseRep';
    
    
    [vIdx, inVIdx] = HPR(fv_small.vertices,poseRep(1,:),1.8);
    
    scatter3(poseRep(1,1), poseRep(1,2), poseRep(1,3), 50, 'b', 'filled');
    scatter3(fv_small.vertices(vIdx,1), fv_small.vertices(vIdx,2), fv_small.vertices(vIdx,3), 10, 'g', 'filled');
    scatter3(fv_small.vertices(inVIdx,1), fv_small.vertices(inVIdx,2), fv_small.vertices(inVIdx,3), 10, 'm', 'filled');
    
    
    scatter3(poseRep(:,1), poseRep(:,2), poseRep(:,3), 50, 'b', 'filled');
    
    
    for pi = 1:size(poseRep, 1)
        fprintf('%d/%d\n', pi, size(poseRep, 1));
        curCamPos = poseRep(pi, :)';
        curCamPos2 = curCamPos;%RMat*curCamPos+TMat;
%         scatter3(curCamPos2(1), curCamPos2(2), curCamPos2(3), 50, 'm', 'filled');
%         curCamPos2 = curCamPos2 + [0.5 0.5 0.5]';
        remainMixtureIdx = find(isObserved~=1);
        for i=remainMixtureIdx
%             i
             EX = X(:,i) + (X(:,i) - curCamPos2).*1000;
             [points, pos, faceInds] = intersectLineMesh3d([curCamPos2' EX'], fv_small.vertices, fv_small.faces, 0);
             if ~isempty(points)
%                  points
%                  curCamPos2
                  dist = bsxfun(@minus, points', curCamPos2);
                  distance =  sqrt(dist(1,:).^2+dist(2,:).^2+dist(3,:).^2);
                  distMixture = curCamPos2-X(:,i);
                  distanceMixture = sqrt(distMixture'*distMixture);
                  [mv, mi] = min(distance);
    %               mv
                  if mv+tol > distanceMixture
    %                    plot3([curCamPos2(1) points(mi,1)]', [curCamPos2(2) points(mi,2)]', [curCamPos2(3) points(mi,3)]', 'r', 'lineWidth', 3);
    %                    scatter3(points(mi,1), points(mi,2), points(mi,3), 50, 'r', 'filled');
    %                    plot3([curCamPos2(1) X(1,i)]', [curCamPos2(2) X(2,i)]', [curCamPos2(3) X(3,i)]', 'b--');  
%                        scatter3(X(1,i), X(2,i), X(3,i), 30, 'g', 'filled');
                      isObserved(i) = 1; %Observed
                  else
    %                   plot3([curCamPos2(1) points(mi,1)]', [curCamPos2(2) points(mi,2)]', [curCamPos2(3) points(mi,3)]', 'g', 'lineWidth', 3);
%                       scatter3(points(mi,1), points(mi,2), points(mi,3), 50, 'r', 'filled');
%                       plot3([points(mi,1) X(1,i)]', [points(mi,2) X(2,i)]', [points(mi,3) X(3,i)]', 'r');  
%                       scatter3(X(1,i), X(2,i), X(3,i), 30, 'r', 'filled');
                      isObserved(i) = 0; %Blocked
                  end


    %                scatter3(X(1,i), X(2,i), X(3,i), 50, 'm', 'filled');
               
             else
                  %plot3([curCamPos2(1) X(1,i)]', [curCamPos2(2) X(2,i)]', [curCamPos2(3) X(3,i)]', 'k');
%                   scatter3(X(1,i), X(2,i), X(3,i), 50, 'k', 'filled');
                 %points
%                  nHit = nHit+1;
             end
             

        end
    end
    ocl(j, :) = isObserved; 
end
% nHit
for i=find(isObserved==1)
     scatter3(X(1,i), X(2,i), X(3,i), 30, 'g', 'filled');
end

for i=find(isObserved==0)
    scatter3(X(1,i), X(2,i), X(3,i), 30, 'r', 'filled');
end
for i=find(isObserved==2)
    scatter3(X(1,i), X(2,i), X(3,i), 30, 'k', 'filled');
end

saveas(h, sprintf('%s/Observable_%d', resultPath, DataSetId));
close(h);

end

fprintf('%Observed :%d/%d\n', sum(isObserved==1), K);
fprintf('%Blocked :%d/%d\n', sum(isObserved==0), K);
fprintf('%Unknown :%d/%d\n', sum(isObserved==2), K);
return;


DataSetId = 11;
figure;
fv_alligned.Faces = fv{DataSetId}.small.faces;
fv_alligned.Vertices = fv{DataSetId}.small.vertices;
patch(fv_alligned,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');
axis equal;
hold on;

poseMat = cell2mat(Pose_Set{DataSetId})';
% poseMat = bsxfun(@plus, RMat*poseMat', TMat)';
scatter3(poseMat(:,1), poseMat(:,2), poseMat(:,3), 50, 'm', 'filled')
[poseRep] = kmeans(poseMat', 30);
poseRep = poseRep';
scatter3(poseRep(:,1), poseRep(:,2), poseRep(:,3), 50, 'b', 'filled');
    
% 
% axis vis3d;
% saveas(h, sprintf('%s/after_reg', resultPath));
% close(h);
% % detect and remove "bad" centers and "unreliable" points 
% [TVrefined,Xrefined,Xrem] = removePointsAndCenters(TV,X,S,a);
% 
% 
% 
% % figure; scatter3(Xrem(1,:)', Xrem(2,:)', Xrem(3,:)');
% % figure;cellfun(@(TVrefined,clrmap,marker,mkSize) scatter3(TVrefined(1,:),TVrefined(2,:),TVrefined(3,:),mkSize,clrmap,marker),TVrefined,clrmap,marker,markerSize, 'UniformOutput', false);
%     
% % visualize TVrefined.
% h=figure;
% hold on, grid on
% 
%     title('Final registration with unreliable points removed','fontweight','bold','fontsize',12);
%     
%     
%     hg3 = cellfun(@(TVrefined,clrmap,marker,mkSize) scatter3(TVrefined(1,:),TVrefined(2,:),TVrefined(3,:),mkSize.*1.5,[0 0 0],'filled'),TVrefined,clrmap,marker,markerSize, 'UniformOutput', false);
%     hg3 = cellfun(@(TVrefined,clrmap,marker,mkSize) scatter3(TVrefined(1,:),TVrefined(2,:),TVrefined(3,:),mkSize,clrmap,'filled'),TVrefined,clrmap,marker,markerSize, 'UniformOutput', false);
%     
%     
%     legend(strIdx{:});
%     
%     % use the same axes as in the registration process
% %     set(gca,'XLim',XLim,'YLim',YLim,'ZLim',Zlim,'fontweight','bold','children', hg3);
%     set(gca,'XLim',XLim,'YLim',YLim,'ZLim',Zlim,'fontweight','bold');
%     
% %     set(3,'position',get(1,'position')+[0 -510 0 0]);
%     
%     view([40 54]) 
% hold off
% saveas(h, sprintf('%s/final_remove_outlier', resultPath));
% close(h);
% 
% % Visualize bad centers (orange) and good centers (blue).
% h=figure;
% hold on, grid on
% 
%     title('Final GMM means.','fontweight','bold','fontsize',12);
%     
%     scatter3(Xrefined(1,:),Xrefined(2,:),Xrefined(3,:),8,[0 .38 .67],'s');
%     
%     scatter3(Xrem(1,:),Xrem(2,:),Xrem(3,:),40,[1 .1412 0],'marker','x');
%     
%     legend('"Good" Centers','"Bad" Centers');
%     
%     % use the same axes as in the registration process
%     set(gca,'XLim',XLim,'YLim',YLim,'ZLim',Zlim,'fontweight','bold');
%     
%     %set(4,'position',get(1,'position')+[+580 -510 0 0]);
%     
%     view([40 54])
%     
% hold off
% saveas(h, sprintf('%s/final_gmm', resultPath));
% close(h);
% 
% 
% clusterSize = 25;
% var_threshold = 1.5;
% [ A, maxIdx, assignments, GMM_mean_color, A_binary, A_binary_c, isStatic ] = posteriorAnalysis( a, X, C2, K, clusterSize, var_threshold );
% 
% 
% id=1;
% gmmk=find(isStatic);
% h=plotGMMK( id, gmmk, assignments, maxIdx, fv, I2, 2);
% 
% %
% h=figure; plot(A'); title(sprintf('Group size of each GMM (#K : %d)', K));
% saveas(h, sprintf('%s/groupassign_k%d', resultPath, K));
% close(h);
% 
% 
% return;
% 
% 
% for id=1:d_size
%     fv3.Faces = fv{id}.Faces;
%     fv3.Vertices = fv{id}.Vertices;
% 
%     
%     sColor = C2{id};
%     mIdx = maxIdx{id};
%     sRMean2 = accumarray(mIdx', sColor(1,:)', [], @mean);
%     sGMean2 = accumarray(mIdx', sColor(2,:)', [], @mean);
%     sBMean2 = accumarray(mIdx', sColor(3,:)', [], @mean);
% 
%     GMM_color = [sRMean2 sGMean2 sBMean2];
%     sampleColor = GMM_color(mIdx, :)';
% 
%     sIndex = I2{id};
%     fv3.FaceVertexCData = sampleColor(:,sIndex)';
% 
%     h = figure;
%     patch(fv3,'FaceColor','flat','EdgeColor','flat',...
%           'MarkerFaceColor','flat', ...
%           'FaceLighting',    'gouraud',     ...
%              'AmbientStrength', 0.15);
%     title('color from GMM');     
%     axis vis3d;   
%     saveas(h, sprintf('%s/color_GMM_%d', resultPath, id));
%     close(h);
% 
%     sampleColor = rGMM_color(mIdx, :)';
%     sIndex = I2{id};
%     fv3.FaceVertexCData = sampleColor(:,sIndex)';
% 
%     h = figure;
%     patch(fv3,'FaceColor','flat','EdgeColor','flat',...
%              'EdgeColor',       'none',        ...
%              'FaceLighting',    'gouraud',     ...
%              'AmbientStrength', 0.15);
%     title('color from GMM');     
%     camlight('headlight');
%     material('dull');
%     axis vis3d;   
%     saveas(h, sprintf('%s/rcolor_GMM_%d', resultPath, id));
%     close(h);
%     
% 
%     fv2.Faces =  fv{id}.Faces;
%     fv2.Vertices = fv{id}.Vertices;
%     sColor = C2{id};
%     fv2.FaceVertexCData = sColor(:,I2{id})';
%     % 
%     h = figure;
%     patch(fv2,'FaceColor','flat','EdgeColor','flat',...
%           'MarkerFaceColor','flat', ...
%           'FaceLighting',    'gouraud',     ...
%              'AmbientStrength', 0.15);
%       title('color from Sampling');     
%     axis vis3d;
%     saveas(h, sprintf('%s/color_sampling_%d.fig', resultPath, id));
%     close(h);
% end
% 
% 
% h = drawKmeansVisTrend( A', assignments, clusterSize );
% % saveas(h, sprintf('%s/groupassign_kmeans%s', resultPath, labels));
% % close(h);
% 
% h = drawKmeansVisTrendBinary( A_binary, assignments, clusterSize );
% saveas(h, sprintf('%s/groupassign_kmeans_binary%s', resultPath, labels));
% close(h);


%Refinement

[R2,t2,X2,S2,a2,pk2,T2,TAssigned2, TXQ2] = jrmpc_soft(V,X,'maxNumIter',refineIterMax, 'updatepriors', 0, 'r', R, 't', t, 's', S, 'initialpriors', pk', 'updateVis', 0);

for epsilonST = [0.05 0.01 0.005 0.001 0.0005 0.0001]
refineIterMax = 50;
% [R2,t2,X2,S2,a2,pk2,T2] = jrmpc_vis(V,Xin,'maxNumIter',refineIterMax,'gamma',0.1, 'updatepriors', 0, 'r', R, 't', t, 's', S, 'vis', A_binary_c);


% [R3,t3,X3,S3,a3,pk3,T3,TAssigned3, TXQ3] = jrmpc_soft(V,X,'maxNumIter',refineIterMax, 'updatepriors', 0, 'r', R, 't', t, 's', S, 'initialpriors', pk', 'updateVis', 0, 'vis', A_binary_c);
% [R4,t4,X4,S4,a4,pk4,T4,TAssigned4, TXQ4] = jrmpc_soft(V,X,'maxNumIter',refineIterMax, 'updatepriors', 0, 'r', R, 't', t, 's', S, 'initialpriors', pk', 'updateTR', 0);
[R5,t5,X5,S5,a5,pk5,T5,TAssigned5, TXQ5] = jrmpc_soft(V,X,'maxNumIter',refineIterMax, 'updatepriors', 0, 'r', R, 't', t, 's', S, 'initialpriors', pk', 'epsilonST', epsilonST);


%Without vis term, just 30 more refine with similar setting

ajk2 = cellfun(@(a) sum(a), a2, 'uniformoutput', false);
ajkMat2 = cell2mat(ajk2);

[ st2 ] = genUniformDist( ajkMat2, epsilonST );
vis2 = genVisFromPeriod(st5, M, K, epsilonST);
vis2 = cell2mat(vis2');

for id=1:2:10
clear 'paramsPlotGMMK';
paramsPlotGMMK.setId = id;
paramsPlotGMMK.gmmk = [];
paramsPlotGMMK.assignments = [];
paramsPlotGMMK.maxIdx = TAssigned2;
paramsPlotGMMK.fv = fv;
paramsPlotGMMK.I2 = I2;
paramsPlotGMMK.type = 2;
paramsPlotGMMK.gmmk_gmmidx = find( (st2(:,2)-st2(:,1)) == M-1 );

h=plotGMMKv2( paramsPlotGMMK );
title(sprintf('Without vis Term, setID=%d, e=%.4f', id, epsilonST));
saveas(h, sprintf('%s/WO_%d_%.4f.png', resultPath, id, epsilonST));
close(h);
end

%With vis term, just 30 more refine with similar setting
% epsilonST = 0.05;
ajk5 = cellfun(@(a) sum(a), a5, 'uniformoutput', false);
ajkMat5 = cell2mat(ajk5);

[ st5 ] = genUniformDist( ajkMat5, epsilonST );
vis5 = genVisFromPeriod(st5, M, K, epsilonST);
vis5 = cell2mat(vis5');

for id=1:2:10
clear 'paramsPlotGMMK';
paramsPlotGMMK.setId = id;
paramsPlotGMMK.gmmk = [];
paramsPlotGMMK.assignments = [];
paramsPlotGMMK.maxIdx = TAssigned5;
paramsPlotGMMK.fv = fv;
paramsPlotGMMK.I2 = I2;
paramsPlotGMMK.type = 2;
paramsPlotGMMK.gmmk_gmmidx = find( (st5(:,2)-st5(:,1)) == M-1 );

h=plotGMMKv2( paramsPlotGMMK );
title(sprintf('With vis Term, setID=%d, e=%.4f', id, epsilonST));
saveas(h, sprintf('%s/W_%d_%.4f.png', resultPath, id, epsilonST));
close(h);
end

end




% 
% 
% params.interval = 1;
% params.marker = marker;
% params.markerSize = markerSize;
% params.clrmap = clrmap;
% params.strIdx = strIdx;
% h = figure;
% drawTransformation(V, T4, params);
% 
% 
% 
% params.type = 2;
% params.Assigned = TAssigned5;
% params.K = K;
% params.interval = 1;
% params.marker = marker;
% params.markerSize = markerSize;
% params.clrmap = clrmap;
% params.strIdx = strIdx;
% 
% h = figure;
% [TV] = drawTransformation(V, T5, params);
% axis equal;
% 
% 
% [ A2, maxIdx2, assignments2, GMM_mean_color2, A_binary2, A_binary_c2, isStatic2 ] = posteriorAnalysis( a2, X2, C2, K, clusterSize, var_threshold );
% 
% 
% epsilon = 0.05;
% ajk = cellfun(@(a) sum(a), a5, 'uniformoutput', false);
% ajkMat = cell2mat(ajk);
% 
% [ st ] = genUniformDist( ajkMat, epsilon );
% vis = genVisFromPeriod(st, M, K);
% vis = cell2mat(vis');
% 
% h = drawKmeansVisTrendBinary( vis, assignments2, clusterSize );
% h = drawKmeansVisTrendBinary( A_binary2, assignments2, clusterSize );
% 
% h = drawKmeansVisTrendBinary( A_binary2-vis, assignments2, clusterSize );
% mean(mean(A_binary2-vis))
% 
% 
% h = drawKmeansVisTrendBinary( vis_old_mat, assignments2, clusterSize );
% h = drawKmeansVisTrendBinary( vis_new_mat-vis_old_mat, assignments2, clusterSize );
% 
% 
% 
% [R5,t5,X5,S5,a5,pk5,T5,TAssigned5, TXQ5] = jrmpc_soft(V,X,'maxNumIter',refineIterMax, 'updatepriors', 0, 'r', R, 't', t, 's', S, 'initialpriors', pk', 'updateVis', 0, 'vis', A_binary_c);

% 
% maxIdx = {};
% 
% for ii=1:length(a2)
%     [m, i] = max(a3{ii}');
%     maxIdx{ii} = i;
%     gmmCount = accumarray(i', 1);
%     A22(1:length(gmmCount), ii) = gmmCount;
% end
% 
% h=figure; plot(A22'); title(sprintf('Group size of each GMM (#K : %d)', K));


% for id=1%:d_size
%     fv3.Faces = fv{id}.Faces;
%     fv3.Vertices = fv{id}.Vertices;
% 
%     
% %     sColor = C2{id};
% %     mIdx = maxIdx{id};
% %     sRMean2 = accumarray(mIdx', sColor(1,:)', [], @mean);
% %     sGMean2 = accumarray(mIdx', sColor(2,:)', [], @mean);
% %     sBMean2 = accumarray(mIdx', sColor(3,:)', [], @mean);
% % 
% %     GMM_color = [sRMean2 sGMean2 sBMean2];
% %     sampleColor = GMM_color(mIdx, :)';
% % 
% %     sIndex = I2{id};
% %     fv3.FaceVertexCData = sampleColor(:,sIndex)';
% % 
% %     h = figure;
% %     patch(fv3,'FaceColor','flat','EdgeColor','flat',...
% %           'MarkerFaceColor','flat', ...
% %           'FaceLighting',    'gouraud',     ...
% %              'AmbientStrength', 0.15);
% %     title('color from GMM');     
% %     axis vis3d;   
% %     saveas(h, sprintf('%s/color_GMM_%d', resultPath, id));
% %     close(h);
% 
%     sampleColor = rGMM_color(mIdx, :)';
%     sIndex = I2{id};
%     fv3.FaceVertexCData = sampleColor(:,sIndex)';
% 
%     h = figure;
%     patch(fv3,'FaceColor','flat','EdgeColor','flat',...
%              'EdgeColor',       'none',        ...
%              'FaceLighting',    'gouraud',     ...
%              'AmbientStrength', 0.15);
%     title('color from GMM');     
%     camlight('headlight');
%     material('dull');
%     axis vis3d;   
% %     saveas(h, sprintf('%s/rcolor_GMM_%d', resultPath, id));
% %     close(h);
%     
% 
% %     fv2.Faces =  fv{id}.Faces;
% %     fv2.Vertices = fv{id}.Vertices;
% %     sColor = C2{id};
% %     fv2.FaceVertexCData = sColor(:,I2{id})';
% %     % 
% %     h = figure;
% %     patch(fv2,'FaceColor','flat','EdgeColor','flat',...
% %           'MarkerFaceColor','flat', ...
% %           'FaceLighting',    'gouraud',     ...
% %              'AmbientStrength', 0.15);
% %       title('color from Sampling');     
% %     axis vis3d;
% %     saveas(h, sprintf('%s/color_sampling_%d.fig', resultPath, id));
% %     close(h);
% end


% [ A2, maxIdx2, assignments2, GMM_mean_color2, A_binary2, A_binary_c2, isStatic2 ] = posteriorAnalysis( a2, X2, C2, K, clusterSize, var_threshold );
% 
% id=1;
% gmmk=find(isStatic2);
% paramsPlotGMMK.setId = id;
% paramsPlotGMMK.gmmk = gmmk;
% paramsPlotGMMK.assignments = assignments2;
% paramsPlotGMMK.maxIdx = TAssigned2;
% paramsPlotGMMK.fv = fv;
% paramsPlotGMMK.I2 = I2;
% paramsPlotGMMK.type = 2;
% paramsPlotGMMK.gmmk_gmmidx = [];
% 
% h=plotGMMKv2( paramsPlotGMMK );
% 
% h = drawKmeansVisTrend( A2', assignments2, clusterSize );
% a = axes;
% t1 = title('K-means based Visability Term');
% a.Visible = 'off'; % set(a,'Visible','off');
% t1.Visible = 'on'; % set(t1,'Visible','on');
% % saveas(h, sprintf('%s/groupassign_kmeans%s.png', resultPath, labels));
% % close(h);
% 
% h = drawKmeansVisTrendBinary( A_binary2, assignments2, clusterSize );
% a = axes;
% t1 = title('K-means based Visability Term');
% a.Visible = 'off'; % set(a,'Visible','off');
% t1.Visible = 'on'; % set(t1,'Visible','on');
% % saveas(h, sprintf('%s/groupassign_kmeans_binary%s.png', resultPath, labels));
% % close(h);
% 
% 
% for epsilon = 0.001%[0.05 0.01 0.005 0.001 0.0005 0.0001]
% ajk2 = cellfun(@(a) sum(a), a2, 'uniformoutput', false);
% ajkMat2 = cell2mat(ajk2);
% [ st ] = genUniformDist( ajkMat2, epsilon );
% vis = genVisFromPeriod(st, 10, K);
% visMat = cell2mat(vis');
% h = drawKmeansVisTrendBinary( visMat, assignments2, clusterSize );  
% a = axes;
% t1 = title(sprintf('Probabilistic Visability Term, e=%f', epsilon));
% a.Visible = 'off'; % set(a,'Visible','off');
% t1.Visible = 'on'; % set(t1,'Visible','on');
% % saveas(h, sprintf('%s/groupassign_prob_binary_e%.4f%s.png', resultPath, epsilon, labels));
% % close(h);
% end
% 
% 
% 
% [ A3, maxIdx3, assignments3, GMM_mean_color3, A_binary3, A_binary_c3, isStatic3 ] = posteriorAnalysis( a3, X3, C2, K, clusterSize, var_threshold );
% gmmk=find(isStatic3);
% 
% id=1;
% h=plotGMMK( id, gmmk, assignments3, TAssigned3, fv, I2, 2);

