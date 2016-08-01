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
resultPath = 'results/EF_vis2';
dataRoot = '/home/mhlee/Kinect_Logs_test';
data_set = {'2016-04-29.00_', '2016-04-29.00_'};
d_size = length(data_set);

[dataSet, poseSet] = loadEFDataset( dataRoot, data_set);
dataSetSampled = loadUniformSampling(dataSet, 20);

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

V = {};
N = {};
C = {};
for i=1:length(dataSet)
    V{i} = dataSet{i}.vs;
    N{i} = dataSet{i}.ns;
    C{i} = dataSet{i}.cs;
end
V = V';
N = N';
C = C';



% return;




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
    markerSize{end+1} = (1+ii)*5;
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

% initialize GMM means Xin, using random sampling of a unit sphere. Choose
% your own initialization. You may want to initialize Xin with some of the
% sets.

% set K as the 50% of the median cardinality of the views
K = ceil(0.5*median(cellfun(@(V) size(V,2),V))); 
K = 50;

rGMM_color = bsxfun(@plus, rand(K, 3)./2, [0.3 0.3 0.3]);

% sample the unit sphere, by randomly selecting azimuth / elevation angles
az = 2*pi*rand(1,K);
el = 2*pi*rand(1,K);


randIdx = randperm(size(V{1}', 1));
% K = 50;
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
set(gca,'fontweight','bold');

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
tic


% [R,t,X,S,a,pk,T,TAssigned, TXQ] = jrmpc_soft_with_normal(V,Xin,'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0);
maxNumIter = 200;
[R,t,X,S,a,pk,T,TAssigned, TXQ] = jrmpc_soft_with_normal(V,Xin, ...
    'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0, ...
    'normal', N, 'normalLambda', .1, 'color', C, 'colorLambda', .1);


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
% axis equal;

return;


[R,t,X,S,a,pk,T,TAssigned, TXQ] = jrmpc_soft_with_normal(V,Xin,'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0);
toc
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




% r_between = cell2mat(T(1,1,end))'*cell2mat(T(2,1,end));
% t_between = -cell2mat(T(1,2,end))+cell2mat(T(2,2,end));
% 
% syn_r_between = syn_r{1}'*syn_r{2}';
% err_r = sqrt(sum(sum(((syn_r_between) - (r_between)).^2)));
% syn_t_between = -syn_t{1}+syn_t{2};
% err_t_1 = sqrt(sum(sum(((syn_t_between)-(t_between)).^2)));
% syn_t_between = syn_t{1}-syn_t{2};
% err_t_2 = sqrt(sum(sum(((syn_t_between)-(t_between)).^2)));
% 
% err_t = min(err_t_1, err_t_2);
% 
% err = err_r + err_t;
% fprintf('err :%f\n', err);

% return;