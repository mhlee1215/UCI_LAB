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
% clear all
% g = gpuDevice(1);
% reset(g);

addpath(genpath('libs'));
data_set = [1 2 3 4];
d_size = length(data_set);

% given an angle theta in radians angle2rotation(theta) construtcts a 
% rotation matrix, rotating a vector with angle theta along all three axes.
% Right hand convention is used for the rotation around the perpendicular.
angle2rotation = @(theta) [1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)] ...
                         *[cos(theta) 0 sin(theta);0 1 0;-sin(theta) 0 cos(theta)] ...
                         *[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];
 
% number of iterations to be run
maxNumIter = 100;                    
                     
% latent angles, rotating V{j} by theta(j) reprojects it to V{1} rotated by
% theta(1). used to quantify the accuracy of the estimated R
theta = [0; pi/20; pi/10; pi/6];
theta = theta(1:d_size);

% number of views, M files must be found in the directory ./syntheticData 
M = numel(theta);

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

sampleNum = 5000;
fv = {};
V2 = {};
V_all = {};


dist_set = {};
dist_set{end+1} = 0.4;
dist_set{end+1} = 0.4;
dist_set{end+1} = 0.4;
dist_set{end+1} = 0.4;

R = [0.36 0.48 -0.8 ; -0.8 0.6 0 ; 0.48 0.64 0.6];

normVar = [];
% minP = [];

for data_idx = data_set%1:d_size
    data_path = sprintf('data/desk/%d', data_idx);
    data_name = 'MeshedReconstruction';
    data_ext = 'stl';
    crop_path = sprintf('%s/crop', data_path);
    fv{data_idx} = stlread(sprintf('%s/%s_%.2f.%s', crop_path, data_name, dist_set{data_idx}, data_ext));
%     rand_idx = randperm(size(fv{data_idx}.vertices, 1));
%     V2{end+1} = fv{data_idx}.vertices(rand_idx(1:sampleNum), :)';

    vertices = fv{data_idx}.vertices;
%     minP = min(vertices);
%     maxP = max(vertices);
%     
%     qFactor = 100;
%     qVertices = uint16((vertices - repmat(minP, size(vertices, 1), 1)).*qFactor)+1;
%     
%     sXMean = accumarray(qVertices, vertices(:,1), [], @mean);
%     sYMean = accumarray(qVertices, vertices(:,2), [], @mean);
%     sZMean = accumarray(qVertices, vertices(:,3), [], @mean);
%     
%     validIdx = find((sXMean ~= 0).*(sYMean ~= 0).*(sZMean ~= 0));
%     vertices_sub = [sXMean(validIdx), sYMean(validIdx), sZMean(validIdx)];
    
%     V2{end+1} = randSubSample(fv{data_idx}.vertices, sampleNum);
    [sData] = uniformSubSample(fv{data_idx}.vertices, 50);
    
%     fvs = reducepatch(fv{data_idx}, 0.01);
    
    if data_idx == 1
        normVar = sqrt(var(sData')');
    end
%     sData = sData ./ repmat(normVar, 1, size(sData, 2));
    if data_idx > 1
%          sData = R*sData;% + repmat(rand(3, 1), 1, size(sData, 2));
%          sData = sData + repmat(rand(3, 1), 1, size(sData, 2));
    end
    V2{end+1} = sData;
%     V2{end+1} = fvs.vertices';
    
    V_all{end+1} = fv{data_idx}.vertices';
end
V = V2';
% return;



% ground truth rotation matrices Rgt{j}*V{j} is aligned with Rgt{1}*R{1}
Rgt = arrayfun(@(theta) angle2rotation(theta),theta,'uniformoutput',false);

% colors for each view
clrmap = {[1 .1412 0]; [.1373 .4196 .5569]; [0 0 1]; [.8039 .6078 .1137]};
clrmap = clrmap(1:d_size);

% markerSizes 
markerSize = {7; 70; 12; 10};
markerSize = markerSize(1:d_size);

% markers
marker = {'s'; 'x'; '.'; '^'};
marker = marker(1:d_size);

% initialize GMM means Xin, using random sampling of a unit sphere. Choose
% your own initialization. You may want to initialize Xin with some of the
% sets.

% set K as the 50% of the median cardinality of the views
K = 300;%ceil(0.5*median(cellfun(@(V) size(V,2),V))); 

% sample the unit sphere, by randomly selecting azimuth / elevation angles
az = 2*pi*rand(1,K);
el = 2*pi*rand(1,K);

%points on a unit sphere
Xin = [cos(az).*cos(el); sin(el); sin(az).*cos(el)];% (unit) polar to cartesian conversion

Xin = Xin/10; % it is good for the initialization to have initial cluster centers at the same order with the points
% since sigma is automatically initialized based on X and V

meanData = cellfun(@(a) mean(a')',V,'uniformoutput',false); % 1 x K rows
meanData = mean(cat(3,meanData{:}),3);
% Xin = Xin ./ repmat(sqrt(var(Xin'))', 1, size(Xin, 2));
Xin = Xin + repmat(meanData, 1, size(Xin, 2));


% show the initial position of the point clouds
figure(1);
hold on,grid on

% make the legend
title('Initial position of the point clouds','fontweight','bold','fontsize',12);

hg1 = cellfun(@(V,clrmap,marker,markerSize) scatter3(V(1,:),V(2,:),V(3,:),markerSize,clrmap,marker),V,clrmap,marker,markerSize, 'UniformOutput', false);

legend(strIdx{:});

set(1,'position',get(1,'position')+[-260 0 0 0]);

% set(gca,'fontweight','bold','children',hg1);
set(gca,'fontweight','bold');

view([40 54])
scatter3(Xin(1,:),Xin(2,:),Xin(3,:),'k')
hold off; drawnow

fprintf('Data registration... \n\n');

tic;
% call JRMPC (type jrmpc with no arguments to see the documentation).
% [R,t,X,S,a,pk,T] = jrmpc(V,Xin,'maxNumIter',maxNumIter,'gamma',0.1);
% [R,t,X,S,a,pk,T] = jrmpc(V,Xin,'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 1);

[R,t,X,S,a,pk,T] = jrmpc2(V,Xin,'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0);
% toc;
% tic;
% [R,t,X,S,a,pk,T] = jrmpc2_gpu(V,Xin,'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 1);
toc;

% measure and display convergency, view 1 is ommited as is the referential.
fprintf('                  ||Rgt{j} - R{j}^T*R{1}||_F                  \n');

fprintf('______________________________________________________________\n');

fprintf('Set  :'),for j=2:M,fprintf('    %d    ',j),end,fprintf('\n');

fprintf('Error:'),for j=2:M,fprintf('  %.4f ',norm(Rgt{j}-R{j}'*R{1},'fro'));end

fprintf('\n');


% visualize the registration process, see documentation of jrmpc for T.

iter=maxNumIter;
TV_all = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V_all',T(:,1,iter),T(:,2,iter),'uniformoutput',false);

figure(2);

for iter = 1:5:maxNumIter
    % apply transformation of iteration : iter
    TV = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V,T(:,1,iter),T(:,2,iter),'uniformoutput',false);
    
    clf(2);
    
    hold on, grid on
    
    title(sprintf('Registration of the sets after %d iteration(s).\n',iter),'fontweight','bold','fontsize',12);
    
    hg2 = cellfun(@(TV,clrmap,marker,markerSize) scatter3(TV(1,:),TV(2,:),TV(3,:),markerSize,clrmap,marker), TV, clrmap, marker, markerSize, 'UniformOutput', false);
    
    legend(strIdx{:});
    
    set(2,'position',get(1,'position')+[+580 0 0 0]);
    
    % iteration 1 locks the axes of subsequent plots
    if iter == 1
       XLim = get(gca,'XLim');
       
       YLim = get(gca,'YLim');
       
       Zlim = get(gca,'ZLim');
       
%        set(gca,'fontweight','bold','children',hg2);
         set(gca,'fontweight','bold');
    else
%        set(gca,'XLim',XLim,'YLim',YLim,'ZLim',Zlim,'fontweight','bold','children',hg2); 
         set(gca,'XLim',XLim,'YLim',YLim,'ZLim',Zlim,'fontweight','bold'); 
    end

    view([40 54])
    
    hold off
    
    pause(.12);
end

% detect and remove "bad" centers and "unreliable" points 
[TVrefined,Xrefined,Xrem] = removePointsAndCenters(TV,X,S,a);



% figure; scatter3(Xrem(1,:)', Xrem(2,:)', Xrem(3,:)');
% figure;cellfun(@(TVrefined,clrmap,marker,mkSize) scatter3(TVrefined(1,:),TVrefined(2,:),TVrefined(3,:),mkSize,clrmap,marker),TVrefined,clrmap,marker,markerSize, 'UniformOutput', false);
    
% visualize TVrefined.
figure(3);
hold on, grid on

    title('Final registration with unreliable points removed','fontweight','bold','fontsize',12);
    
    hg3 = cellfun(@(TVrefined,clrmap,marker,mkSize) scatter3(TVrefined(1,:),TVrefined(2,:),TVrefined(3,:),mkSize,clrmap,marker),TVrefined,clrmap,marker,markerSize, 'UniformOutput', false);
    
    legend(strIdx{:});
    
    % use the same axes as in the registration process
%     set(gca,'XLim',XLim,'YLim',YLim,'ZLim',Zlim,'fontweight','bold','children', hg3);
    set(gca,'XLim',XLim,'YLim',YLim,'ZLim',Zlim,'fontweight','bold');
    
    set(3,'position',get(1,'position')+[0 -510 0 0]);
    
    view([40 54]) 
hold off

% Visualize bad centers (orange) and good centers (blue).
figure(4);
hold on, grid on

    title('Final GMM means.','fontweight','bold','fontsize',12);
    
    scatter3(Xrefined(1,:),Xrefined(2,:),Xrefined(3,:),8,[0 .38 .67],'s');
    
    scatter3(Xrem(1,:),Xrem(2,:),Xrem(3,:),40,[1 .1412 0],'marker','x');
    
    legend('"Good" Centers','"Bad" Centers');
    
    % use the same axes as in the registration process
    set(gca,'XLim',XLim,'YLim',YLim,'ZLim',Zlim,'fontweight','bold');
    
    set(4,'position',get(1,'position')+[+580 -510 0 0]);
    
    view([40 54])
    
hold off

sqe = @(Y,X) sum(bsxfun(@minus,permute(Y,[2 3 1]),permute(X,[3 2 1])).^2,3);
ddd = cellfun(@(TV) sqe(TV,X),TV,'uniformoutput',false);
ddd2 = mean([ddd{1} ; ddd{2}],1);
ddd2 = mean([ddd{1}],1);

figure;scatter3(X(1,:),X(2,:),X(3,:),40,ddd2,'filled'); colorbar;
figure;scatter3(X(1,:),X(2,:),X(3,:),40,pk,'filled'); colorbar;
figure;scatter3(X(1,:),X(2,:),X(3,:),40,S,'filled'); colorbar;


A = [];
A2 = [];
A3 = [];
A4 = [];%zeros(K, 1);
for ii=1:d_size
[m, i] = max(a{ii}');
A = [A accumarray(i', 1)];
A2 = [A2 sum(a{ii})'];
A3 = [A3 mean(a{ii})'];

[sortedValues,sortIndex] = sort(a{ii}','descend'); 


firstIdx = sortIndex(1, :)';
firstVal = sortedValues(1, :)';
secondIdx = sortIndex(2, :)';
secondVal = sortedValues(2, :)';

KKK = accumarray(firstIdx, firstVal);
KKK2 = accumarray(secondIdx, secondVal);

A4 = [A4 KKK+KKK2];
end
AA = bsxfun(@rdivide, A, max(A')');
% figure; plot(AA');
% figure; plot(A');
% figure; plot(A2');
% figure; plot(A3');
% figure; plot(A4');

gmm_score = AA;

cluster = [];
for i=1:K
    cluster = [cluster kmeans(gmm_score(i,:)', 2)];
end
figure; plot(cluster)

