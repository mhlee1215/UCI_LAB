function [ bestR, bestT ] = featureBasedRegistration( v_m1, c_m1, n_m1, v_m2, c_m2, n_m2, d_seeds_1, f_1, d_seeds_2, f_2, trial, isConstrained )
%FEATUREBASEDREGISTRATION Summary of this function goes here
%   Detailed explanation goes here
% 

if nargin < 12
    isConstrained = false;
end


if ~isempty(c_m1)
    [v_m1_s,~,c_m1_s] = uniformSubSample(v_m1, 15, c_m1);
    v_m1_s = v_m1_s';
    c_m1_s = c_m1_s';
else
    c_m1_s = c_m1;
end

if ~isempty(c_m2)
    [v_m2_s,~,c_m2_s] = uniformSubSample(v_m2, 15, c_m2);
    v_m2_s = v_m2_s';
    c_m2_s = c_m2_s';
else
    c_m2_s = c_m2;
end

if ~isempty(n_m1)
    [v_m1_s,~,n_m1_s] = uniformSubSample(v_m1, 15, n_m1);
    v_m1_s = v_m1_s';
    n_m1_s = n_m1_s';
else
    n_m1_s = n_m1;
end

if ~isempty(n_m2)
    [v_m2_s,~,n_m2_s] = uniformSubSample(v_m2, 15, n_m2);
    v_m2_s = v_m2_s';
    n_m2_s = n_m2_s';
else
    n_m2_s = n_m2;
end



% v_m1_s = v_m1;
% c_m1_s = c_m1;
% n_m1_s = n_m1;
% v_m2_s = v_m2;
% c_m2_s = c_m2;
% n_m2_s = n_m2;

if ~isempty(c_m1_s)
    c_m1_s = gaussianNormalize(c_m1_s);
end
% if ~isempty(n_m1_s)
%     n_m1_s = gaussianNormalize(n_m1_s);
% end

if ~isempty(c_m2_s)
    c_m2_s = gaussianNormalize(c_m2_s);
end
% if ~isempty(n_m2_s)
%     n_m2_s = gaussianNormalize(n_m2_s);
% end



% figure; hold on; 
% scatter3(v_m1(:,1), v_m1(:,2), v_m1(:,3), 'k', 'filled'); 
% scatter3(d_seeds_1(1,:), d_seeds_1(2,:), d_seeds_1(3,:), 'r', 'filled');



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


% view_matches( vertex11, c_m1_s, vertex22, c_m2_s, key11, key22, NN_Data, aa, 0.99);
params.compAll = 1;
params.density = 30; %higher 
[dist, seeds, overlap, d] = voxelizer_compare(v_m1_s, c_m1_s, n_m1_s, v_m2_s, c_m2_s, n_m2_s, params); 
[seeds_valid_1] = getValidSeeds(v_m1_s, seeds, d);
[seeds_valid_2] = getValidSeeds(v_m2_s, seeds, d);
minOverlap = min(size(seeds_valid_1, 1), size(seeds_valid_2, 1))*0.3;

validRatio = 0.95;
validIdx = find((aa(:,1) ./ aa(:,2) < validRatio));
disp(sprintf('valid length :%d', length(validIdx)));
% trial = 300;
minDist = dist
bestR = [1 0 0 ; 0 1 0 ; 0 0 1];
bestT = [0 0 0]';

if length(validIdx) < 3
    bestR = [];
    bestT = [];
    return;
end

for t=1:trial
    if mod(t, 10) == 0
        fprintf('trial .. %d/%d\n', t, trial);
    end
    %sample matches random way
    matchNum = 3;
    topK = 1;
    randIdx = randperm(length(validIdx));
    randValidIdx = validIdx(randIdx(1:matchNum));
    A = d_seeds_1(:, randValidIdx);
    
    [~, randI] = max(rand(topK, matchNum));
    b = [];
    for pp = 1:matchNum
        b = [b d_seeds_2(:, NN_Data(randValidIdx(pp),1))];
    end
    %Get rigid transform
    [R, T] = rigid_transform_3D(b', A');
    
    if isConstrained == true
        [ R, T ] = fcn_supressConstraint( R, T );
    end
%     R_degree = sum(sum((R-eye(3)).^2));
%     if R_degree > 0.8
%         continue;
%     end
    
    v_m2_s_t = bsxfun(@plus, R*v_m2_s', T)';
%     tic

    
%     n_m1_s = (Rn1*n_m1_s')';
    n_m2_s_t = (R*n_m2_s')';

%     figure; hold on; 
%     scatter3(v_m1_s(:,1), v_m1_s(:,2), v_m1_s(:,3), 'k', 'filled'); 
%     scatter3(v_m2_s_t(:,1), v_m2_s_t(:,2), v_m2_s_t(:,3), 'r', 'filled'); 
%     
%     figure; hold on; 
%     scatter3(v_m1_s(:,1), v_m1_s(:,2), v_m1_s(:,3), 'k', 'filled'); 
% %     scatter3(v_m2_s(:,1), v_m2_s(:,2), v_m2_s(:,3), 'r', 'filled'); 
%     scatter3(d_seeds_1(1,:), d_seeds_1(2,:), d_seeds_1(3,:), 'r', 'filled');

    
    [dist, seeds, overlap, d] = voxelizer_compare(v_m1_s, c_m1_s, n_m1_s, v_m2_s_t, c_m2_s, n_m2_s_t, params); 
%     [dist2, seeds2, overlap2, d2] = voxelizer_compare(v_m1_s, c_m1_s, n_m1_s, v_m2_s, c_m2_s, n_m2_s, params); 
%     pclviewer([v_m1_s c_m1_s ; v_m2_s_t c_m2_s]');
    %dist = dist% / overlap;
%     overlap
%     toc
    
%     sprintf('%f, %f', minDist, dist);
    if minDist > dist
        minDist = dist;
        bestR = R;
        bestT = T;
        fprintf('minDist: %f\n', dist);
        fprintf('bestR:\n'); bestR
        fprintf('bestT:\n'); bestT
    end
end

% [dist, seeds, overlap, d] = voxelizer_compare(v_m1_s, c_m1_s, v_m2_s_best, c_m2_s, params); 

% v_m2_best = bsxfun(@plus, bestR*v_m2', bestT)';
% pclviewer([v_m1 c_m1 ; v_m2_best c_m2]');


end

