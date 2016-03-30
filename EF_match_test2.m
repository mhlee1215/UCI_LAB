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


v_m1 = fv{1}.Vertices;
v_m2 = fv{2}.Vertices;
c_m1 = fv{1}.FaceVertexCData;
c_m2 = fv{2}.FaceVertexCData;

[v_m1_s,~,c_m1_s] = uniformSubSample(v_m1, 5, c_m1);
v_m1_s = v_m1_s';
c_m1_s = c_m1_s';
[v_m2_s,~,c_m2_s] = uniformSubSample(v_m2, 5, c_m2);
v_m2_s = v_m2_s';
c_m2_s = c_m2_s';


%Subsampling for dense descriptors
params_desc.normalRadius=0.2;
params_desc.searchRadius=1.2;
params_desc.searchK = 0;

[d_seeds_1,~,~] = uniformSubSample(v_m1_s, 5, c_m1_s);
f_1 = FPFHSDescriptor(v_m1_s, c_m1_s, d_seeds_1, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);
mean(sum(f_1==0,1))
 
[d_seeds_2,~,~] = uniformSubSample(v_m2_s, 5, c_m2_s);
f_2 = FPFHSDescriptor(v_m2_s, c_m2_s, d_seeds_2, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);


[ bestR, bestT ] = featureBasedRegistration( v_m1_s, c_m1_s, v_m2_s, c_m2_s, d_seeds_1, f_1, d_seeds_2, f_2, 300 );
v_m2_best11 = bsxfun(@plus, bestR*v_m2_s', bestT)';
pclviewer([v_m1_s c_m1_s ; v_m2_best11 c_m2_s]');

%%After registration without warping
v_m2_best22 = bsxfun(@plus, bestR*v_m2', bestT)';
pclviewer([v_m1 c_m1 ; v_m2_best22 c_m2]');

%%Before registration without warping
pclviewer([v_m1 c_m1 ; v_m2 c_m2]');

pclviewer([v_m1 c_m1]');
pclviewer([v_m2 c_m2]');


%% Do part registration
%For key 1
key1 = SIFTKeypoint(v_m1_s, c_m1_s, '');
%Remove duplicated key
key1 = unique(key1', 'rows')';

%For key 2
key2 = SIFTKeypoint(v_m2_best11, c_m2_s, '');
%Remove duplicated key
key2 = unique(key2', 'rows')';
f_1_key = FPFHSDescriptor(v_m1_s, c_m1_s, key1, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);
f_2_key = FPFHSDescriptor(v_m2_best11, c_m2_s, key2, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);


DescrData = pdist2(f_1_key',f_2_key');
posDistance = pdist2(key1', key2');
searchRange = 0.8;
posFilteredDescData = DescrData + (posDistance > searchRange).*999;

[pos_sorted_dist,pos_NN_Data]=sort(posFilteredDescData,2,'ascend');
pos_NN_Data=pos_NN_Data(:,1:NN_Number);





%% Clustering
validIdx = find((pos_sorted_dist(:,1) ./ pos_sorted_dist(:,2) < 0.7) ...
                .* (pos_sorted_dist(:,1) ./ pos_sorted_dist(:,2) > 0.1));

gKey1 = double(key1(:, validIdx));
fv1_v = fv{1}.Vertices;
fv1_c = fv{1}.FaceVertexCData;
IDX_1 = knnsearch(gKey1', fv1_v);

%% Soft clustering
kdtree = vl_kdtreebuild(gKey1) ;
[index, distance] = vl_kdtreequery(kdtree, gKey1, fv1_v', 'NumNeighbors', 5, 'MaxComparisons', 55) ;


% pclviewer([gKey1']');

gKey2 = key2(:, pos_NN_Data(validIdx, 1));
fv2_v = fv{2}.Vertices;
fv2_v = bsxfun(@plus, bestR*fv2_v', bestT)';
fv2_c = fv{2}.FaceVertexCData;
IDX_2 = knnsearch(gKey2', fv2_v);






%% Match between clusters
% ii_set = [6, 8, 11, 13, 15, 16, 17];
bestSubRSet = {};
bestSubTSet = {};

IDX_SET = unique(IDX_1);
IDX_SET_invalid = [];
for iidx = 1:length(IDX_SET)
    ii = IDX_SET(iidx);
    ii_idx = find(IDX_1==ii);
    ii_idx2 = find(IDX_2==ii);

    if isempty(ii_idx) || isempty(ii_idx2)
        IDX_SET_invalid = [ IDX_SET_invalid ; ii];
        continue;
    end

    v_m1 = [fv1_v(ii_idx, :)];
    c_m1 = [fv1_c(ii_idx, :)];
    % pclviewer([v_m1 c_m1]');

    v_m2 = [fv2_v(ii_idx2, :)];
    c_m2 = [fv2_c(ii_idx2, :)];
    % pclviewer([v_m2 c_m2]');


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

    [ bestSubR, bestSubT ] = featureBasedRegistration( v_m1_s2, c_m1_s2, v_m2_s2, c_m2_s2, d_seeds_12, f_12, d_seeds_22, f_22, 400 );

    % v_m2_best = bsxfun(@plus, bestSubR*v_m2', bestSubT)';
    % pclviewer([v_m1 c_m1 ; v_m2_best c_m2]');
    % pclviewer([v_m1 c_m1 ; v_m2 c_m2]');


    bestSubRSet{ii} = bestSubR;
    bestSubTSet{ii} = bestSubT;

end


v_m2_best = bsxfun(@plus, bestSubR*v_m2', bestSubT)';
pclviewer([v_m1 c_m1 ; v_m2_best c_m2]');
pclviewer([v_m1 c_m1 ; v_m2 c_m2]');

% v_m2_s_best = bsxfun(@plus, bestSubR*v_m2_s', bestSubT)';
% 
% pclviewer([v_m1_s ;v_m2_s_best]');

%% View Result
iidx = 10;
ii = IDX_SET(iidx);
ii_idx = find(IDX_1==ii);
ii_idx2 = find(IDX_2==ii);

v_m1 = [fv1_v(ii_idx, :)];
c_m1 = [fv1_c(ii_idx, :)];
% pclviewer([v_m1 c_m1]');

v_m2 = [fv2_v(ii_idx2, :)];
c_m2 = [fv2_c(ii_idx2, :)];
% pclviewer([v_m2 c_m2]');
v_m2_best = bsxfun(@plus, bestSubRSet{ii}*v_m2', bestSubTSet{ii})';
pclviewer([v_m1 c_m1 ; v_m2_best c_m2]');

pclviewer([v_m1 c_m1 ; v_m2 c_m2]');

pclviewer([v_m1 c_m1]');
pclviewer([v_m2 c_m2]');


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
viewMergedModel(fv1_v, fv1_c, fv2_v, fv2_c, bestSubRSet, bestSubTSet, IDX_SET, IDX_1);



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


