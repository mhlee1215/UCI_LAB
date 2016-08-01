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


%data_set = [1:6 8:14];%[1 2 3 4];
dataRoot = '/home/mhlee/data';
% data_set = {'2016-02-09.01','2016-02-09.02','2016-02-09.03'};
data_set = {'2016-02-24.00','2016-02-24.02'};
d_size = length(data_set);

resultPath = 'results/EF_portable_1';
mkdir(resultPath);


maxNumIter = 100;                    
                     

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


NN_Number = 5;

% sampleNum = 5000;
fv = {};
V2 = {};
C2 = {};
I2 = {};
V_all = {};
Pose_Set = {};
PoseMat_Set = {};


dist_set = {};
dist_set{end+1} = 0.4;
dist_set{end+1} = 0.4;
dist_set{end+1} = 0.4;
dist_set{end+1} = 0.4;

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

return;

v_m1_all = fv{1}.Vertices;
v_m2_all = fv{2}.Vertices;

c_m1_all = fv{1}.FaceVertexCData;
c_m2_all = fv{2}.FaceVertexCData;

n_m1_all = fv{1}.Normal;
n_m2_all = fv{2}.Normal;

[v_m1_s,~,c_m1_s] = uniformSubSample(v_m1_all, 5, c_m1_all);
v_m1_s = v_m1_s';
c_m1_s = c_m1_s';
[v_m1_s,~,n_m1_s] = uniformSubSample(v_m1_all, 5, n_m1_all);
v_m1_s = v_m1_s';
n_m1_s = n_m1_s';

[v_m2_s,~,c_m2_s] = uniformSubSample(v_m2_all, 5, c_m2_all);
v_m2_s = v_m2_s';
c_m2_s = c_m2_s';
[v_m2_s,~,n_m2_s] = uniformSubSample(v_m2_all, 5, n_m2_all);
v_m2_s = v_m2_s';
n_m2_s = n_m2_s';



%Subsampling for dense descriptors
params_desc.normalRadius=0.2;
params_desc.searchRadius=1.2;
params_desc.searchK = 0;

[d_seeds_1,~,~] = uniformSubSample(v_m1_s, 5, c_m1_s);
f_1 = FPFHSDescriptor(v_m1_s, c_m1_s, d_seeds_1, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);
mean(sum(f_1==0,1))
 
[d_seeds_2,~,~] = uniformSubSample(v_m2_s, 5, c_m2_s);
f_2 = FPFHSDescriptor(v_m2_s, c_m2_s, d_seeds_2, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);


[ initBestR, initBestT ] = featureBasedRegistration( v_m1_s, c_m1_s, n_m1_s, v_m2_s, c_m2_s, n_m2_s, d_seeds_1, f_1, d_seeds_2, f_2, 300 );
v_m2_best11 = bsxfun(@plus, initBestR*v_m2_s', initBestT)';
pclviewer([v_m1_s c_m1_s ; v_m2_best11 c_m2_s]');

%%After registration without warping
% [ bestR, bestT ] = featureBasedRegistration( v_m1_s, c_m1_s, [], v_m2_s, c_m2_s, [], d_seeds_1, f_1, d_seeds_2, f_2, 300 );
% [ bestR, bestT ] = featureBasedRegistration( v_m1_s, c_m1_s, n_m1_s, v_m2_s, c_m2_s, n_m2_s, d_seeds_1, f_1, d_seeds_2, f_2, 300 );
% v_m2_best22 = bsxfun(@plus, bestR*v_m2_all', bestT)';
% pclviewer([v_m1_all c_m1_all ; v_m2_best22 c_m2_all]');

%%Before registration without warping
pclviewer([v_m1_all c_m1_all ; v_m2_all c_m2_all]');

pclviewer([v_m1_all c_m1_all]');
pclviewer([v_m2_all c_m2_all]');


%% Do part registration


%Dense sampling but has less radious for
samplingDensity = 5;
[v_m1_s,~,c_m1_s] = uniformSubSample(v_m1_all, samplingDensity, c_m1_all);
v_m1_s = v_m1_s';
c_m1_s = c_m1_s';
[v_m1_s,~,n_m1_s] = uniformSubSample(v_m1_all, samplingDensity, n_m1_all);
v_m1_s = v_m1_s';
n_m1_s = n_m1_s';

[v_m2_s,~,c_m2_s] = uniformSubSample(v_m2_all, samplingDensity, c_m2_all);
v_m2_s = v_m2_s';
c_m2_s = c_m2_s';
[v_m2_s,~,n_m2_s] = uniformSubSample(v_m2_all, samplingDensity, n_m2_all);
v_m2_s = v_m2_s';
n_m2_s = n_m2_s';



%For key 1
key1 = SIFTKeypoint(v_m1_s, c_m1_s, '');
%Remove duplicated key
key1 = unique(key1', 'rows')';

%For key 2
key2 = SIFTKeypoint(v_m2_best11, c_m2_s, '');
%Remove duplicated key
key2 = unique(key2', 'rows')';

samplingDensity = 30;
[v_m1_s,~,c_m1_s] = uniformSubSample(v_m1_all, samplingDensity, c_m1_all);
v_m1_s = v_m1_s';
c_m1_s = c_m1_s';
[v_m1_s,~,n_m1_s] = uniformSubSample(v_m1_all, samplingDensity, n_m1_all);
v_m1_s = v_m1_s';
n_m1_s = n_m1_s';

[v_m2_s,~,c_m2_s] = uniformSubSample(v_m2_all, samplingDensity, c_m2_all);
v_m2_s = v_m2_s';
c_m2_s = c_m2_s';
[v_m2_s,~,n_m2_s] = uniformSubSample(v_m2_all, samplingDensity, n_m2_all);
v_m2_s = v_m2_s';
n_m2_s = n_m2_s';

params_desc.normalRadius=0.2;
params_desc.searchRadius=0.2;
params_desc.searchK = 0;

v_m2_best11 = bsxfun(@plus, initBestR*v_m2_s', initBestT)';

f_1_key = FPFHSDescriptor(v_m1_s, c_m1_s, key1, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);
f_2_key = FPFHSDescriptor(v_m2_best11, c_m2_s, key2, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);


DescrData = pdist2(f_1_key',f_2_key');
posDistance = pdist2(key1', key2');
searchRange = 0.8;
posFilteredDescData = DescrData + (posDistance > searchRange).*999;

[pos_sorted_dist,pos_NN_Data]=sort(posFilteredDescData,2,'ascend');
pos_NN_Data=pos_NN_Data(:,1:NN_Number);





%% Clustering
fv1_v = fv{1}.Vertices;
fv1_c = fv{1}.FaceVertexCData;
fv1_n = fv{1}.Normal;

fv2_v = fv{2}.Vertices;
fv2_v = bsxfun(@plus, initBestR*fv2_v', initBestT)';
fv2_c = fv{2}.FaceVertexCData;
fv2_n = fv{2}.Normal;


validIdx = find((pos_sorted_dist(:,1) ./ pos_sorted_dist(:,2) < 0.7) ...
                .* (pos_sorted_dist(:,1) ./ pos_sorted_dist(:,2) > 0.1));

            
% [d_seeds_all,~,~] = uniformSubSample([fv1_v ;fv2_v] , 4, [fv1_c ; fv2_c]);
            
gKey1 = double(key1(:, validIdx));
gKey2 = double(key2(:, pos_NN_Data(validIdx, 1)));


group = {};
group{1} = Clustering(fv1_v, gKey1, 'soft');
group{2} = Clustering(fv2_v, gKey2, 'soft');

% return;
%% sub clustering

clustering_id = 2;

fv_sub = {};
clustersCells = {};
clustersCenters = {};
w = [1.0 0.0 0.0];
group_soft = {};

for g_id = 1:length(fv)
    g_id
    fv_sub{g_id}.v = fv{g_id}.Vertices(group{g_id}{clustering_id}, :);
    fv_sub{g_id}.c = fv{g_id}.FaceVertexCData(group{g_id}{clustering_id}, :);
    fv_sub{g_id}.n = fv{g_id}.Normal(group{g_id}{clustering_id}, :);
    % pclviewer([fv_1_v_s fv_1_c_s]');
    
%     fv_v_s_n = gaussianNormalize(fv_sub{g_id}.v);
%     fv_c_s_n = gaussianNormalize(fv_sub{g_id}.c);
%     fv_n_s_n = gaussianNormalize(fv_sub{g_id}.n);
    
    fv_v_s_n = fv_sub{g_id}.v;%gaussianNormalize(fv_sub{g_id}.v);
    fv_c_s_n = fv_sub{g_id}.c;%gaussianNormalize(fv_sub{g_id}.c);
    fv_n_s_n = fv_sub{g_id}.n;%gaussianNormalize(fv_sub{g_id}.n);
    
   
    fv_v_s_n(find(isnan(fv_v_s_n))) = 0;
    fv_c_s_n(find(isnan(fv_c_s_n))) = 0;
    fv_n_s_n(find(isnan(fv_n_s_n))) = 0;

    [clustCent,point2cluster,clustMembsCell] = ...
        MeanShiftCluster([fv_v_s_n.*w(1) fv_c_s_n.*w(2) fv_n_s_n.*w(3)]', 0.2);
    %cluster from meanshift
    clustersCells{g_id} = clustMembsCell;
    %centers
    clustersCenters{g_id} = clustCent(1:3, :);
    %soft clustering
    group_soft{g_id} = Clustering(fv_v_s_n, clustersCenters{g_id}, 'soft');
end

g_id = 1;
figure; hold on;
for sub_g_id = 1:length(clustersCells{g_id})
g_idx = clustersCells{g_id}{sub_g_id};%find(point2cluster==g_id);
scatter3(fv_sub{g_id}.v(g_idx, 1), fv_sub{g_id}.v(g_idx, 2), fv_sub{g_id}.v(g_idx, 3), 1, 'filled');
end

g_id_src = 1;
g_id_dst = 2;

C_src = clustersCenters{g_id_src};
C_dst = clustersCenters{g_id_dst};



indexKNN = getKNN(C_src, C_dst);

RSet = {};
TSet = {};

params.normalRadius=0.05;
params.searchRadius=0.1;
params.searchK = 0;
params.dataUnitSize = 15;
params.keypointUnitSize = 15;
params.descriptorUnitSize = 100;
params.d_seeds_12 = [];
params.d_seeds_22 = [];
params.f_12 = [];
params.f_22 = [];


for idx_cluster_src = 5%2:size(C_src, 2)
    
    ID_src = clustersCells{g_id_src}{idx_cluster_src};
        
    fv_sub{g_id_src}.v = fv{g_id_src}.Vertices(group{g_id_src}{clustering_id}, :);
    fv_sub{g_id_src}.c = fv{g_id_src}.FaceVertexCData(group{g_id_src}{clustering_id}, :);
    fv_sub{g_id_src}.n = fv{g_id_src}.Normal(group{g_id_src}{clustering_id}, :);

    v1 = fv_sub{g_id_src}.v(ID_src, :);
    c1 = fv_sub{g_id_src}.c(ID_src, :);
    n1 = fv_sub{g_id_src}.n(ID_src, :);
    
    v1 = fv_sub{g_id_src}.v;
    c1 = fv_sub{g_id_src}.c;
    n1 = fv_sub{g_id_src}.n;
    
    
%     for idx_cluster_dst2 = 1:length(indexKNN(:, idx_cluster_src))
%         idx_cluster_dst = indexKNN(idx_cluster_dst2, idx_cluster_src);
%         if idx_cluster_src == 0 || idx_cluster_dst == 0
%             continue;
%         end
%         
%         ID_dst = group_soft{g_id_dst}{idx_cluster_dst};
%         
%         fv_sub{g_id_dst}.v = fv{g_id_dst}.Vertices(group{g_id_dst}{clustering_id}, :);
%         fv_sub{g_id_dst}.c = fv{g_id_dst}.FaceVertexCData(group{g_id_dst}{clustering_id}, :);
%         fv_sub{g_id_dst}.n = fv{g_id_dst}.Normal(group{g_id_dst}{clustering_id}, :);
% 
%         v2 = fv_sub{g_id_dst}.v(ID_dst, :);
%         c2 = fv_sub{g_id_dst}.c(ID_dst, :);
%         n2 = fv_sub{g_id_dst}.n(ID_dst, :);
        
        v2 = fv_sub{g_id_dst}.v;
        c2 = fv_sub{g_id_dst}.c;
        n2 = fv_sub{g_id_dst}.n;
        
        
        
        
        
        [R, T, d_seeds_12, f_12, d_seeds_22, f_22 ] = featureBasedRegistrationWraper(v1, c1, n1, v2, c2, n2, params);
        
        
        params.d_seeds_22 = d_seeds_22;
        params.f_22 = f_22;
        
%         RSet{idx_cluster_src}{idx_cluster_dst2} = R;
%         TSet{idx_cluster_src}{idx_cluster_dst2} = T;
        
        RSet{idx_cluster_src} = R;
        TSet{idx_cluster_src} = T;
%         break;
%         pclviewer([v1 c1 ; v2 c2]');
%         v1_best = (R'*bsxfun(@minus, v1', T))';
%         pclviewer([v1_best c1 ; v2 c2]');
        
%     end
%     break;
end






fv2 = [];
fv2.Vertices = fv_sub{g_id}.v(1:1000,:);
tri = delaunayTriangulation(fv2.Vertices(:,1), fv2.Vertices(:,2), fv2.Vertices(:,3));
tri2 = delaunay(fv2.Vertices(:,1), fv2.Vertices(:,2), fv2.Vertices(:,3));
fv2.Faces = tri2(:,2:4);
fv2.FaceVertexCData = [1 0 0];

% figure;
% scatter3(fv2.Vertices(:,1), fv2.Vertices(:,2), fv2.Vertices(:,3), 0.2);

patch(fv2, ...
             'FaceColor','flat','EdgeColor','flat', ...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);
        camlight('headlight');
        material('dull');
         axis vis3d;


%  patch(fv2,'FaceColor','flat','EdgeColor','flat',...
%   'MarkerFaceColor','flat', ...
%   'FaceLighting',    'gouraud',     ...
%      'AmbientStrength', 0.15);
% 
% 
% 




g_id_src = 1;
g_id_dst = 2;
idx_cluster_src = 3;

ID_src = clustersCells{g_id_src}{idx_cluster_src};
        
fv_sub{g_id_src}.v = fv{g_id_src}.Vertices(group{g_id_src}{clustering_id}, :);
fv_sub{g_id_src}.c = fv{g_id_src}.FaceVertexCData(group{g_id_src}{clustering_id}, :);
fv_sub{g_id_src}.n = fv{g_id_src}.Normal(group{g_id_src}{clustering_id}, :);

v1 = fv_sub{g_id_src}.v(ID_src, :);
c1 = fv_sub{g_id_src}.c(ID_src, :);
n1 = fv_sub{g_id_src}.n(ID_src, :);

v2 = fv_sub{g_id_dst}.v;
c2 = fv_sub{g_id_dst}.c;
n2 = fv_sub{g_id_dst}.n;

R = RSet{idx_cluster_src};
T = TSet{idx_cluster_src};
v1_best = (R'*bsxfun(@minus, v1', T))';

pclviewer([v1_best c1 ; v2 c2]');
% pclviewer([v1 c1]');
% pclviewer([v2 c2]');


g_id = 2;
figure; hold on;
idx_dst = unique(index_filtered);
for sub_g_id2 = 2:length(idx_dst)%length(clustersCells{g_id})
    sub_g_id = idx_dst(sub_g_id2);
    g_idx = clustersCells{g_id}{sub_g_id};%find(point2cluster==g_id);
scatter3(fv_sub{g_id}.v(g_idx, 1), fv_sub{g_id}.v(g_idx, 2), fv_sub{g_id}.v(g_idx, 3), 1, 'filled');
end

%Data ID
g_id1 = 1;
g_id2 = 2;

%Sub Group ID
sub_g_id_1 = 1;
sub_g_id_2 = 1;

g_idx = {};
g_idx{g_id1} = clustersCells{g_id1}{sub_g_id_1};
g_idx{g_id2} = clustersCells{g_id2}{sub_g_id_2};

v_m1 = fv_sub{g_id1}.v(g_idx{g_id1}, :);
c_m1 = fv_sub{g_id1}.c(g_idx{g_id1}, :);
n_m1 = fv_sub{g_id1}.n(g_idx{g_id1}, :);

v_m2 = fv_sub{g_id2}.v(g_idx{g_id2}, :);
c_m2 = fv_sub{g_id2}.c(g_idx{g_id2}, :);
n_m2 = fv_sub{g_id2}.n(g_idx{g_id2}, :);

pclviewer([v_m1 c_m1 ; v_m2 c_m2]');

pclviewer([v_m1 c_m1]');
pclviewer([v_m2 c_m2]');

[bestSubR, bestSubT ] = featureBasedRegistrationWraper(v_m1, c_m1, n_m1, v_m2, c_m2, n_m2);
v_m2_best = bsxfun(@plus, bestSubR*v_m2', bestSubT)';
pclviewer([v_m1 c_m1 ; v_m2_best c_m2]');

pclviewer([v_m1 c_m1]');
pclviewer([v_m1 c_m1 ; v_m2 c_m2]');




%% Match between clusters
% ii_set = [6, 8, 11, 13, 15, 16, 17];
bestSubRSet = {};
bestSubTSet = {};

%IDX_SET = unique(IDX_1);
IDX_SET_invalid = [];
for ii = 19%1:length(group{1})%length(IDX_SET)
    %ii = IDX_SET(iidx);
    
    %get group
    ii_idx = group{1}{ii};%find(IDX_1==ii);
    ii_idx2 = group{2}{ii};%find(IDX_2==ii);

    if isempty(ii_idx) || isempty(ii_idx2)
        IDX_SET_invalid = [ IDX_SET_invalid ; ii];
        continue;
    end

    v_m1 = fv{1}.Vertices(ii_idx, :);
    c_m1 = fv{1}.FaceVertexCData(ii_idx, :);
%     % pclviewer([v_m1 c_m1]');
% 
    v_m2 = fv{2}.Vertices(ii_idx2, :);
    v_m2 = bsxfun(@plus, initBestR*v_m2', initBestT)';
    c_m2 = fv{2}.FaceVertexCData(ii_idx2, :);
%     

    [v_m1_s2,~,c_m1_s2] = uniformSubSample(v_m1, 15, c_m1);
    v_m1_s2 = v_m1_s2';
    c_m1_s2 = c_m1_s2';
    [v_m2_s2,~,c_m2_s2] = uniformSubSample(v_m2, 15, c_m2);
    v_m2_s2 = v_m2_s2';
    c_m2_s2 = c_m2_s2';


    %Subsampling for dense descriptors
    params_desc.normalRadius=0.2;
    params_desc.searchRadius=0.8;
    params_desc.searchK = 0;

    [d_seeds_12,~,~] = uniformSubSample(v_m1_s2, 4, c_m1_s2);
    f_12 = FPFHSDescriptor(v_m1_s2, c_m1_s2, d_seeds_12, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);
    % mean(sum(f_12==0, 1))
    [d_seeds_22,~,~] = uniformSubSample(v_m2_s2, 4, c_m2_s2);
    f_22 = FPFHSDescriptor(v_m2_s2, c_m2_s2, d_seeds_22, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);
    % mean(sum(f_22==0, 1))

    [ bestSubR, bestSubT ] = featureBasedRegistration( v_m1_s2, c_m1_s2, [], v_m2_s2, c_m2_s2, [], d_seeds_12, f_12, d_seeds_22, f_22, 400 );

    
%     v_m2_best = bsxfun(@plus, bestSubR*v_m2', bestSubT)';
%     pclviewer([v_m1 c_m1 ; v_m2_best c_m2]');
% %     pclviewer([v_m1 c_m1 ; v_m2 c_m2]');


    bestSubRSet{ii} = bestSubR;
    bestSubTSet{ii} = bestSubT;

end


% v_m2_best = bsxfun(@plus, bestSubR*v_m2', bestSubT)';
% pclviewer([v_m1 c_m1 ; v_m2_best c_m2]');
% pclviewer([v_m1 c_m1 ; v_m2 c_m2]');

% v_m2_s_best = bsxfun(@plus, bestSubR*v_m2_s', bestSubT)';
% 
% pclviewer([v_m1_s ;v_m2_s_best]');

%% View Result

ii = 9;

%Get group
ii_idx = group_1{ii};
ii_idx2 = group_2{ii};

v_m1 = [fv1_v(ii_idx, :)];
c_m1 = [fv1_c(ii_idx, :)];
% pclviewer([v_m1 c_m1]');

v_m2 = [fv2_v(ii_idx2, :)];
c_m2 = [fv2_c(ii_idx2, :)];
% pclviewer([v_m2 c_m2]');
v_m2_best = bsxfun(@plus, bestSubRSet{ii}*v_m2', bestSubTSet{ii})';
pclviewer([v_m1 c_m1 ; v_m2_best c_m2]');
% 
% pclviewer([v_m1 c_m1 ; v_m2 c_m2]');

pclviewer([v_m1 c_m1]');
% pclviewer(gKey1)
% pclviewer([v_m2 c_m2]');


%%Merge into a single big model Type 1
v_m1_all = fv1_v;
c_m1_all = fv1_c;

v_m2_all = [];
c_m2_all = [];
for iidx=1:length(IDX_SET)
    ii = IDX_SET(iidx);
    
    if ii > length(bestSubRSet) || isempty(bestSubRSet{ii})
        RR = eye(3);
        TT = [0 ; 0 ; 0];
%         continue;
    else
        RR = bestSubRSet{ii};
        TT = bestSubTSet{ii};
    end
    
    
%     ii_idx = find(IDX_1==ii);
    ii_idx2 = find(IDX_2==ii);

%     v_m1 = [fv1_v(ii_idx, :)];
%     c_m1 = [fv1_c(ii_idx, :)];
    % pclviewer([v_m1 c_m1]');

    v_m2 = [fv2_v(ii_idx2, :)];
    c_m2 = [fv2_c(ii_idx2, :)];
    % pclviewer([v_m2 c_m2]');
    v_m2_best = bsxfun(@plus, RR*v_m2', TT)';
    
    v_m2_all = [v_m2_all ; v_m2_best];
    c_m2_all = [c_m2_all ; c_m2];
    
end
pclviewer([v_m2_all c_m2_all]');
pclviewer([v_m1_all c_m1_all;v_m2_all c_m2_all]');


%%Merge into a single big model Type 2
viewMergedModel(fv1_v, fv1_c, fv2_v, fv2_c, bestSubRSet, bestSubTSet, group{1}, group{2});



% Trans = GMBasedRegistration(v_m1_s, c_m1_s, v_m2_s, c_m2_s, 0.0);
% RR = Trans(:,1,end);
% TT = Trans(:,2,end);
% v_m2_best = (RR{1}'*bsxfun(@minus, bsxfun(@plus, RR{2}*v_m2_best', TT{2}), TT{1}))';
% pclviewer([v_m1 c_m1 ; v_m2_best c_m2]');


IDX_standard = IDX_1; %move from 1 to 2
V_standrad = fv1_v;

sXMean = accumarray(IDX_standard, V_standrad(:,1), [], @mean);
sYMean = accumarray(IDX_standard, V_standrad(:,2), [], @mean);
sZMean = accumarray(IDX_standard, V_standrad(:,3), [], @mean);

V_centers = [sXMean(:), sYMean(:), sZMean(:)]';
V_centers(find(V_centers==0)) = -999; %Exclude empty centers
V_centers(:, IDX_SET_invalid) = -999;

kdtree = vl_kdtreebuild(V_centers) ;
[index, distance] = vl_kdtreequery(kdtree, V_centers, V_centers, 'NumNeighbors', 10, 'MaxComparisons', 55) ;

validIdx = distance(find(distance < 100));

bestSubRSet_Q = [];
bestSubTSet_Q = [];
for i=1:length(bestSubRSet)
    if isempty(bestSubRSet{i}) continue; end
    bestSubRSet_Q(i, :) = q_getFromRotationMatrix(bestSubRSet{i});
    bestSubTSet_Q(i, :) = bestSubTSet{i};
end

IDX_standard_set = unique(IDX_standard);


lambda = 0.5;
nn_smooth = 3;

cur_bestSubRSet_Q = bestSubRSet_Q;
cur_bestSubTSet_Q = bestSubTSet_Q;

for i=1:10
    next_bestSubRSet_Q = [];
    next_bestSubTSet_Q = [];
    for i = 1:length(IDX_standard_set)
        IDX_next = IDX_standard_set(i);

        w = [lambda (1-lambda).*(1./(distance(2:nn_smooth,IDX_next)'))];
        subRSet_q_cur = cur_bestSubRSet_Q(index(1:nn_smooth, IDX_next), :);
        ave_r_q = q_getWeightedAverage(subRSet_q_cur, w);
        ave_t_q = sum(repmat(w', 1, 3).*cur_bestSubTSet_Q(index(1:nn_smooth, IDX_next), :), 1);

        next_bestSubRSet_Q(IDX_next, :) = ave_r_q;
        next_bestSubTSet_Q(IDX_next, :) = ave_t_q';
    end
    cur_bestSubRSet_Q = next_bestSubRSet_Q;
    cur_bestSubTSet_Q = next_bestSubTSet_Q;
end


next_bestSubRSet = {};
next_bestSubTSet = {};
for i=1:size(next_bestSubRSet_Q, 1)
    next_bestSubRSet{i} = q_getRotationMatrix(next_bestSubRSet_Q(i,:));
    next_bestSubTSet{i} = next_bestSubTSet_Q(i,:)';
end

viewMergedModel(fv1_v, fv1_c, fv2_v, fv2_c, next_bestSubRSet, next_bestSubTSet, IDX_SET, IDX_1);



rMat = RR{1};
qs = q_getFromRotationMatrix(rMat);
rMat2 = q_getRotationMatrix(qs);


return;

Trans = GMBasedRegistration(v_m1_s, c_m1_s, v_m2_s_best, c_m2_s, 0.0);
RR = Trans(:,1,end);
TT = Trans(:,2,end);
v_m2_best = (RR{1}'*bsxfun(@minus, bsxfun(@plus, RR{2}*v_m2_best', TT{2}), TT{1}))';
pclviewer([v_m1 c_m1 ; v_m2_best c_m2]');

pclviewer([v_m1 c_m1]');
pclviewer([v_m2_best c_m2]');



sum(sum((bestSubRSet{1}-eye(3)).^2))

% save('match_exp_1.mat', 'IDX_SET', 'IDX_1', 'IDX_2', 'fv1_v', 'fv1_c', 'fv2_v', 'fv2_c', 'bestSubRSet', 'bestSubTSet');



gKey1 = key1(:, validIdx);
fv1_v = fv{1}.Vertices;
fv1_c = fv{1}.FaceVertexCData;

IDX_1 = knnsearch(gKey1', fv1_v);


