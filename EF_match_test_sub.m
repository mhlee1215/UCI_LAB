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

load('subData.mat');

[v1_s,~,c1_s] = uniformSubSample(v1, 100, c1);
v1_s = v1_s';
c1_s = c1_s';
[v2_s,~,c2_s] = uniformSubSample(v2, 100, c2);
v2_s = v2_s';
c2_s = c2_s';


params.normalRadius=0.1;
params.searchRadius=0.1;
params.searchK = 0;
[desc1, key1] = getFeatures(v1_s, c1_s, params);
[desc2, key2] = getFeatures(v2_s, c2_s, params);
mean(sum(desc1==0,2))
mean(sum(desc2==0,2))


% FinalPose = icp_mod_point_plane_pyr(v1, n1, v2, n2, 0.02, 200, 3, 1, 8, 0, 1);



vertex1 = v1_s;
color1 = c1_s;
vertex2 = v2_s;
color2 = c2_s;

NN_Number = 5;
DescrData = pdist2(desc1',desc2');
[aa,NN_Data]=sort(DescrData,2,'ascend');
NN_Data=NN_Data(:,1:NN_Number);
% validIdx = find((NN_Data(:,1) ./ NN_Data(:,2) < 0.95));
%[matches, scores] = vl_ubcmatch(desc1, desc2) ;


dist_between = mean(max(v1_s) - min(v1_s) )*2;
vertex11 = vertex1;
vertex11(:,3) = vertex11(:,3).*1;
vertex22 = vertex2+repmat([dist_between 0 0], size(vertex2, 1), 1);
vertex22(:,3) = vertex22(:,3).*1;

key11 = key1;
key11(:,3) = key11(:,3).*1;
key22 = key2+repmat([dist_between 0 0]', 1, size(key2, 2));
key22(:,3) = key22(:,3).*1;

% return;
% 
% view_matches( vertex11, color1, vertex22, color2, key11, key22, NN_Data, aa );
% title('before the first registration');

h=figure;
interval = 1;
for idx = 1:size(key11,2)
    clf;
    hold on;
    view_match_sim( idx, interval, DescrData, vertex11, color1, vertex22, color2, key11, key22, NN_Data, params.searchRadius )
    view(28, -20);
    axis equal;
    pause();

end

return;


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


v_m = [fv_middle{1}.Vertices];
c_m = [fv_middle{1}.FaceVertexCData ];
pclviewer([v_m c_m]');

v_m = [warped];
c_m = [fv_middle{2}.FaceVertexCData];
pclviewer([v_m c_m]');





view_matches( vertex11, color1, vertex22, color2, key11, key22, pos_NN_Data, pos_sorted_dist );
title('after the first registration');


%
validIdx = find((pos_sorted_dist(:,1) ./ pos_sorted_dist(:,2) < 0.7));

gKey1 = key1(:, validIdx);
fv1_v = fv{1}.Vertices;
fv1_c = fv{1}.FaceVertexCData;
IDX_1 = knnsearch(gKey1', fv1_v);

pclviewer([gKey1']');

gKey2 = key2(:, pos_NN_Data(validIdx, 1));
fv2_v = fv{2}.Vertices;
fv2_c = fv{2}.FaceVertexCData;
IDX_2 = knnsearch(gKey2', fv2_v);


ii = 6;
ii_idx = find(IDX_1==ii);
ii_idx2 = find(IDX_2==ii);

v_m1 = [fv1_v(ii_idx, :)];
c_m1 = [fv1_c(ii_idx, :)];
% pclviewer([v_m1 c_m1]');

v_m2 = [fv2_v(ii_idx2, :)];
c_m2 = [fv2_c(ii_idx2, :)];
% pclviewer([v_m2 c_m2]');



featureBasedRegistration(v_m1, c_m1, v_m2, c_m2);

%
%Downsampling
[v_m1_s,~,c_m1_s] = uniformSubSample(v_m1, 15, c_m1);
v_m1_s = v_m1_s';
c_m1_s = c_m1_s';
[v_m2_s,~,c_m2_s] = uniformSubSample(v_m2, 15, c_m2);
v_m2_s = v_m2_s';
c_m2_s = c_m2_s';


%Subsampling for dense descriptors
params_desc.normalRadius=0.2;
params_desc.searchRadius=0.8;
params_desc.searchK = 0;


[d_seeds_1,~,~] = uniformSubSample(v_m1_s, 5, c_m1_s);
f_1 = FPFHSDescriptor(v_m1_s, c_m1_s, d_seeds_1, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);

[d_seeds_2,~,~] = uniformSubSample(v_m2_s, 5, c_m2_s);
f_2 = FPFHSDescriptor(v_m2_s, c_m2_s, d_seeds_2, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);

NN_Number = 5;
DescrData = pdist2(f_1',f_2');
[aa,NN_Data]=sort(DescrData,2,'ascend');
NN_Data=NN_Data(:,1:NN_Number);

dist_between = mean(max(v_m1_s) - min(v_m1_s) )*2;
vertex11 = v_m1_s;
vertex11(:,3) = vertex11(:,3).*1;
vertex22 = v_m2_s+repmat([dist_between 0 0], size(v_m2_s, 1), 1);
vertex22(:,3) = vertex22(:,3).*1;

key11 = d_seeds_1;
key11(:,3) = key11(:,3).*1;
key22 = d_seeds_2+repmat([dist_between 0 0]', 1, size(d_seeds_2, 2));
key22(:,3) = key22(:,3).*1;


view_matches( vertex11, c_m1_s, vertex22, c_m2_s, key11, key22, NN_Data, aa, 0.99);

% params.density = 30; %higher 
% [dist, seeds, overlap, d] = voxelizer_compare(v_m1_s, c_m1_s, v_m2_s, c_m2_s, params); 
% [seeds_valid_1] = getValidSeeds(v_m1_s, seeds, d);
% [seeds_valid_2] = getValidSeeds(v_m2_s, seeds, d);
% minOverlap = min(size(seeds_valid_1, 1), size(seeds_valid_2, 1))*0.3;
% 
% validRatio = 0.99;
% validIdx = find((aa(:,1) ./ aa(:,2) < validRatio));
% trial = 100;
% minDist = 999999;
% bestR = [];
% bestT = [];
% for t=1:trial
%     if mod(t, 10) == 0
%         fprintf('trial .. %d/%d\n', t, trial);
%     end
%     %sample matches random way
%     matchNum = 3;
%     topK = 1;
%     randIdx = randperm(length(validIdx));
%     randValidIdx = validIdx(randIdx(1:matchNum));
%     A = d_seeds_1(:, randValidIdx);
%     
%     [~, randI] = max(rand(topK, matchNum));
%     b = [];
%     for pp = 1:matchNum
%         b = [b d_seeds_2(:, NN_Data(randValidIdx(pp),1))];
%     end
%     %Get rigid transform
%     [R, T] = rigid_transform_3D(b', A');
%     v_m2_s_t = bsxfun(@plus, R*v_m2_s', T)';
% %     tic
%     
%     [dist, seeds, overlap, d] = voxelizer_compare(v_m1_s, c_m1_s, v_m2_s_t, c_m2_s, params); 
%     %dist = dist% / overlap;
% %     toc
%     if minDist > dist && overlap > minOverlap
%         minDist = dist;
%         bestR = R;
%         bestT = T;
%         fprintf('minDist: %f\n', dist);
%     end
% end

% [dist, seeds, overlap, d] = voxelizer_compare(v_m1_s, c_m1_s, v_m2_s_best, c_m2_s, params); 

% 
[v_m1_s,~,c_m1_s] = uniformSubSample(v_m1, 15, c_m1);
v_m1_s = v_m1_s';
c_m1_s = c_m1_s';
[v_m2_s,~,c_m2_s] = uniformSubSample(v_m2, 15, c_m2);
v_m2_s = v_m2_s';
c_m2_s = c_m2_s';


%Subsampling for dense descriptors
params_desc.normalRadius=0.2;
params_desc.searchRadius=0.8;
params_desc.searchK = 0;

[d_seeds_1,~,~] = uniformSubSample(v_m1_s, 5, c_m1_s);
f_1 = FPFHSDescriptor(v_m1_s, c_m1_s, d_seeds_1, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);

[d_seeds_2,~,~] = uniformSubSample(v_m2_s, 5, c_m2_s);
f_2 = FPFHSDescriptor(v_m2_s, c_m2_s, d_seeds_2, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);


[ bestR, bestT ] = featureBasedRegistration( v_m1_s, c_m1_s, v_m2_s, c_m2_s, d_seeds_1, f_1, d_seeds_2, f_2 );
v_m2_best = bsxfun(@plus, bestR*v_m2', bestT)';
pclviewer([v_m1 c_m1 ; v_m2_best c_m2]');

v_m2_s_best = bsxfun(@plus, bestR*v_m2_s', bestT)';

pclviewer([v_m1_s ;v_m2_s_best]');


Trans = GMBasedRegistration(v_m1_s, c_m1_s, v_m2_s_best, c_m2_s, 0.0);
RR = Trans(:,1,end);
TT = Trans(:,2,end);
v_m2_best = (RR{1}'*bsxfun(@minus, bsxfun(@plus, RR{2}*v_m2_best', TT{2}), TT{1}))';
pclviewer([v_m1 c_m1 ; v_m2_best c_m2]');

pclviewer([v_m1 c_m1]');
pclviewer([v_m2_best c_m2]');







