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
data_set = {'2016-04-29.02','2016-04-29.08'};
d_size = length(data_set);

resultPath = 'results/EF_normal_1';
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

fprintf(' Data loading...\n');

% load the views, file view<j>.txt corresponds to theta(j)
% cutView*.txt is a partial view as described in the paper, while view*.txt
% "sees" the whole surface (again downsampled and noisy)

%V = arrayfun(@(j) dlmread(sprintf('./syntheticData/view%d.txt',j),' ')',idx,'uniformoutput',false);
% V = arrayfun(@(j) dlmread(sprintf('libs/JRMPC_v0.9.4/syntheticData/cutView%d.txt',j),' ')',idx,'uniformoutput',false);

% sampleNum = 5000;
fv = {};
V2 = {};
C2 = {};
I2 = {};
V_all = {};
Pose_Set = {};
PoseMat_Set = {};

% R = [0.36 0.48 -0.8 ; -0.8 0.6 0 ; 0.48 0.64 0.6];

normVar = [];
% minP = [];


for data_idx = 1:d_size
%     data_path = sprintf('data/desk/%d', data_idx);
%     data_name = 'MeshedReconstruction';
%     data_ext = 'stl';
%     crop_path = sprintf('%s/crop', data_path);
    
    fprintf('loading.. %d/%d', data_idx, 1);
    
%     fv{data_idx} = stlread(sprintf('%s/%s_%.2f.%s', crop_path, data_name, dist_set{data_idx}, data_ext));
    
    
    data3DFilePath = sprintf('%s/%s.klg.ply', dataRoot, data_set{data_idx});
    [tri, pts, data, comments] = ply_read(data3DFilePath, 'tri');
    
    fprintf('end\n');
    
    v = [data.vertex.x data.vertex.y data.vertex.z];
    f = tri';
    c = [data.vertex.red data.vertex.green data.vertex.blue];
    n = [data.vertex.nx data.vertex.ny data.vertex.nz];
    
%     r_trans = (rand(3,1)-0.5)/2;
%     r_theta = rand()*pi/8 - pi/16;
%     r_rot = [cos(r_theta) -sin(r_theta) 0 ; sin(r_theta) cos(r_theta) 0 ; 0 0 1];
%     
    r_trans = [0 0 0]';
    r_theta = 0;
    r_rot = [cos(r_theta) -sin(r_theta) 0 ; sin(r_theta) cos(r_theta) 0 ; 0 0 1];
    
    
    %Make intentional error
    v = bsxfun(@plus, r_rot*v', r_trans)';
    

    fv{data_idx}.Faces = f;
    fv{data_idx}.Vertices = v;
    fv{data_idx}.FaceVertexCData = c./255;
    fv{data_idx}.Normal = n;
    


    %Create small version to accerelate occlusion finding
%     fprintf('Generate small version..');
%     fv_small.vertices = fv{data_idx}.Vertices;
%     fv_small.faces = fv{data_idx}.Faces;
%     fv_small = reducepatch(fv_small, 0.01);
%     fprintf('End\n');
%     fv{data_idx}.small = fv_small;
    
    
    %Subsampling
    [sData, sIndex, sColor] = uniformSubSample(fv{data_idx}.Vertices, 10, fv{data_idx}.FaceVertexCData);
    
    V2{end+1} = sData;
    C2{end+1} = sColor;
    I2{end+1} = sIndex;
%     Pose_Set{end+1} = cellfun(@(pose) r_rot*pose'+r_trans, poseSet,'uniformoutput',false);
%     PoseMat_Set{end+1} = cellfun(@(pose) [[r_rot*pose(1:3,1:3) ; [0 0 0]] [r_rot*pose(1:end-1,end)+r_trans ;1]], poseMatSet2,'uniformoutput',false);
    V_all{end+1} = fv{data_idx}.Vertices';
end
V = V2';

% return;

fv_small = {};
for data_idx=1:2
    [sData, sIndex, sColor, sNormal] = uniformSubSample(fv{data_idx}.Vertices, 5, fv{data_idx}.FaceVertexCData, fv{data_idx}.Normal);
    fv_small{data_idx}.Vertices = sData';
    fv_small{data_idx}.FaceVertexCData = sColor';
    fv_small{data_idx}.Normal = sNormal';
end

fv_middle = {};
for data_idx=1:2
    [sData, sIndex, sColor, sNormal] = uniformSubSample(fv{data_idx}.Vertices, 20, fv{data_idx}.FaceVertexCData, fv{data_idx}.Normal);
    fv_middle{data_idx}.Vertices = sData';
    fv_middle{data_idx}.FaceVertexCData = sColor';
    fv_middle{data_idx}.Normal = sNormal';
end

fv_target = fv_middle;




% fv_sub = {};
clustersCells = {};
% clustersCenters = {};
% w = [2.0 0.0 1.0];
% group_soft = {};

for g_id = 1:length(fv_target)
    g_id
%     fv_v_s_n = gaussianNormalize(fv_sub{g_id}.v);
%     fv_c_s_n = gaussianNormalize(fv_sub{g_id}.c);
%     fv_n_s_n = gaussianNormalize(fv_sub{g_id}.n);
    
    fv_v_s_n = fv_target{g_id}.Vertices;%gaussianNormalize(fv_sub{g_id}.v);
    fv_c_s_n = fv_target{g_id}.FaceVertexCData;%gaussianNormalize(fv_sub{g_id}.c);
    fv_n_s_n = fv_target{g_id}.Normal;%gaussianNormalize(fv_sub{g_id}.n);
    
   
    fv_v_s_n(find(isnan(fv_v_s_n))) = 0;
    fv_c_s_n(find(isnan(fv_c_s_n))) = 0;
    fv_n_s_n(find(isnan(fv_n_s_n))) = 0;
    
    
    [d_seeds,~,~] = uniformSubSample(fv_v_s_n, 1);
    clustersCells{g_id} = Clustering(fv_v_s_n, d_seeds, 'hard');

%     [clustCent,point2cluster,clustMembsCell] = ...
%         MeanShiftCluster([fv_v_s_n.*w(1) fv_c_s_n.*w(2) fv_n_s_n.*w(3)]', 0.4);
%     %cluster from meanshift
%     clustersCells{g_id} = clustMembsCell;
%     %centers
%     clustersCenters{g_id} = clustCent(1:3, :);
    %soft clustering
%     group_soft{g_id} = Clustering(fv_v_s_n, clustersCenters{g_id}, 'soft');
end


filtered.idx = {};
filtered.threshold = 0.12;
for g_id = 1%:length(fv_target)

    idx_merge = [];
    cvar = [];
    for c_id = 1:length(clustersCells{g_id})
        idx = clustersCells{g_id}{c_id};
        v = fv_target{g_id}.Vertices(idx, :);
        c = fv_target{g_id}.FaceVertexCData(idx, :);
        n = fv_target{g_id}.Normal(idx, :);
        
%         n = find(~isnan(n))

        curVar = mean(var(n));
        cvar = [cvar curVar];
        if curVar < filtered.threshold
            idx_merge = [idx_merge idx];
        end
    end
    

%     figure; histogram(cvar, 10);
        
    filtered.idx{g_id} = idx_merge;

end

randIdx = randperm(size(fv_target{1}.Vertices, 1));
K = 1000;
Xin = fv_target{1}.Vertices(randIdx(1:K),:)';
maxNumIter = 100;                    


tic
[R,t,X,S,a,pk,T,TAssigned, TXQ] = jrmpc_soft(V,Xin,'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0);
toc


drawMarkers;

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

pclviewer([TV{1} TV{2} ;C2{1} C2{2}]);


[clust1, colorP1] = Segmentation(fv{1}.Vertices, fv{1}.FaceVertexCData, 1000, []);
[clust2, colorP2] = Segmentation(fv{2}.Vertices, fv{2}.FaceVertexCData, 1000, []);

k = SIFTKeypoint(fv{2}.Vertices, fv{2}.FaceVertexCData, '0.03 6 10 0.5');

pclviewer([fv{1}.Vertices fv{1}.FaceVertexCData]');


k = SIFTKeypoint(clust1{222}(1:3, :)', clust1{222}(4:6, :)', '0.05 6 10 0.5');




return;











params.normalRadius=0.2;
params.searchRadius=2;
params.searchK = 0;
[desc1, key1] = getFeatures(fv_target{1}.Vertices, fv_target{1}.FaceVertexCData, params);
[desc2, key2] = getFeatures(fv_target{2}.Vertices, fv_target{2}.FaceVertexCData, params);
mean(sum(desc1==0,2))
mean(sum(desc2==0,2))

% return;

vertex1 = fv_target{1}.Vertices;
color1 = fv_target{1}.FaceVertexCData;
vertex2 = fv_target{2}.Vertices;
color2 = fv_target{2}.FaceVertexCData;

NN_Number = 5;
DescrData = pdist2(desc1',desc2');
[aa,NN_Data]=sort(DescrData,2,'ascend');
NN_Data=NN_Data(:,1:NN_Number);

%[matches, scores] = vl_ubcmatch(desc1, desc2) ;


dist_between = mean(max(fv_target{1}.Vertices) - min(fv_target{1}.Vertices) )*2;
vertex11 = vertex1;
vertex11(:,3) = vertex11(:,3).*1;
vertex22 = vertex2+repmat([dist_between 0 0], size(vertex2, 1), 1);
vertex22(:,3) = vertex22(:,3).*1;

key11 = key1;
key11(:,3) = key11(:,3).*1;
key22 = key2+repmat([dist_between 0 0]', 1, size(key2, 2));
key22(:,3) = key22(:,3).*1;


view_matches( vertex11, color1, vertex22, color2, key11, key22, NN_Data, aa );
title('before the first registration');

h=figure;
interval = 1;
for idx = 1:50:size(key11,2)
    clf;
    hold on;
    view_match_sim( idx, interval, DescrData, vertex11, color1, vertex22, color2, key11, key22, NN_Data, params.searchRadius )
    view(28, -20);
    axis equal;
    pause();

end



A = key1';
b = key2(:, NN_Data(:,1))';
[R, t] = rigid_transform_3D(b, A);


fv_1 = fv_target{1};
fv_2 = fv_target{2};

fv_2.Vertices = bsxfun(@plus, R*fv_2.Vertices', t)';
key2_t = bsxfun(@plus, R*key2, t);

%Before
v_m = [fv_1.Vertices ; fv_target{2}.Vertices];
c_m = [fv_1.FaceVertexCData ; fv_2.FaceVertexCData];
pclviewer([v_m c_m]');

%View
v_m = [fv_1.Vertices ; fv_2.Vertices];
c_m = [fv_1.FaceVertexCData ; fv_2.FaceVertexCData];
pclviewer([v_m c_m]');


posDistance = pdist2(key1', key2_t');
searchRange = 0.8;
posFilteredDescData = DescrData + (posDistance > searchRange).*999;

[pos_sorted_dist,pos_NN_Data]=sort(posFilteredDescData,2,'ascend');
pos_NN_Data=pos_NN_Data(:,1:NN_Number);

view_matches( vertex11, color1, vertex22, color2, key11, key22, pos_NN_Data, pos_sorted_dist );
title('after the first registration');

validIdx = find((pos_sorted_dist(:,1) ./ pos_sorted_dist(:,2) < 0.95));

h=figure;
interval = 1;
for idx = 1:10:size(key11,2)
    clf;
    hold on;
    view_match_sim( idx, interval, posFilteredDescData, vertex11, color1, vertex22, color2, key11, key22, pos_NN_Data, searchRange )
    view(28, -20);
    axis equal;
    pause();
end

%Do Warping from evertyhing else to dataset 1

target = key1';
src = key2(:, pos_NN_Data(:,1))';
src_t = bsxfun(@plus, R*src', t)';
ambient = bsxfun(@plus, R*fv_middle{2}.Vertices', t)';
tic;
warped = warp3Dtps2(src_t(validIdx, :), target(validIdx, :),ambient,100);
toc;

figure(1); clf;
plot3(target(validIdx,1),target(validIdx,2),target(validIdx,3),'g.'); hold on;
plot3(src_t(validIdx,1),src_t(validIdx,2),src_t(validIdx,3),'r.');
for ii = 1:length(validIdx)%size(target,1);
    i = validIdx(ii);
  plot3([target(i,1) src_t(i,1)],[target(i,2) src_t(i,2)],[target(i,3) src_t(i,3)],'-');
end

% bsxfun(@plus, R*fv_small{2}.Vertices', t)';
v_m = [fv_middle{1}.Vertices ; ambient];
c_m = [fv_middle{1}.FaceVertexCData ; fv_middle{2}.FaceVertexCData];
pclviewer([v_m c_m]');

v_m = [fv_middle{1}.Vertices ; warped];
c_m = [fv_middle{1}.FaceVertexCData ; fv_middle{2}.FaceVertexCData];
pclviewer([v_m c_m]');


v_m = [warped];
c_m = [fv_middle{2}.FaceVertexCData];
pclviewer([v_m c_m]');










