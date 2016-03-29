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
data_set = 1:8;%[1 2 3 4];
d_size = length(data_set);

resultPath = 'results/r4_3dvar';
mkdir(resultPath);




% return;
     
% ptCloud = pcread('e:\\MeshedReconstruction.ply');
% pcshow(ptCloud);
% 
% 
% [v, f, n, c, stltitle] = stlreadColor('data/desk_color/1/MeshedReconstruction_256.stl');
% [v, f, n, c, stltitle] = stlreadColor('e:\\MeshedReconstruction2.stl');
% figure;
% patch('Vertices', v, 'Faces', f, 'FaceVertexCData', c,...
%          'EdgeColor',       'none',        ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);

% given an angle theta in radians angle2rotation(theta) construtcts a 
% rotation matrix, rotating a vector with angle theta along all three axes.
% Right hand convention is used for the rotation around the perpendicular.
% angle2rotation = @(theta) [1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)] ...
%                          *[cos(theta) 0 sin(theta);0 1 0;-sin(theta) 0 cos(theta)] ...
%                          *[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];
 
% number of iterations to be run
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


dist_set = {};
dist_set{end+1} = 0.4;
dist_set{end+1} = 0.4;
dist_set{end+1} = 0.4;
dist_set{end+1} = 0.4;

% R = [0.36 0.48 -0.8 ; -0.8 0.6 0 ; 0.48 0.64 0.6];

normVar = [];
% minP = [];


for data_idx = data_set%1:d_size
    data_path = sprintf('data/desk/%d', data_idx);
    data_name = 'MeshedReconstruction';
    data_ext = 'stl';
    crop_path = sprintf('%s/crop', data_path);
    
    fprintf('loading.. %d/%d', data_idx, 1);
    
%     fv{data_idx} = stlread(sprintf('%s/%s_%.2f.%s', crop_path, data_name, dist_set{data_idx}, data_ext));
    
    dataRoot = 'E:\lab\work\2015_4d_seg\codes\data\desk_color2';
    [tri, pts, data, comments] = ply_read(sprintf('%s/MeshedReconstruction_%d.ply', dataRoot, data_idx), 'tri');
    
    fprintf('end\n');
    
    
    
    
    v = [data.vertex.x data.vertex.y data.vertex.z];
    f = tri';
    c = [data.vertex.red data.vertex.green data.vertex.blue];

    fv{data_idx}.Faces = f;
    fv{data_idx}.Vertices = v;
    fv{data_idx}.FaceVertexCData = c./255;

    h = figure;
    patch(fv{data_idx},'FaceColor','flat','EdgeColor','flat',...
          'MarkerFaceColor','flat', ...
          'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);
      title('color from Original');    
      axis vis3d;   
    saveas(h, sprintf('%s/color_original_%d.fig', resultPath, data_idx));
    close(h);
%     axis vis3d;
    
    
%     rand_idx = randperm(size(fv{data_idx}.vertices, 1));
%     V2{end+1} = fv{data_idx}.vertices(rand_idx(1:sampleNum), :)';

%     vertices = fv{data_idx}.Vertices;
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
    [sData, sIndex, sColor] = uniformSubSample(fv{data_idx}.Vertices, 70, fv{data_idx}.FaceVertexCData);
    
    
%     fv2.Faces =  fv{data_idx}.Faces;
%     fv2.Vertices = fv{data_idx}.Vertices;
%     sColor = C2{end};
%     fv2.FaceVertexCData = sColor(:,I2{end})';
% % 
%     figure;
%     patch(fv2,'FaceColor','flat','EdgeColor','flat',...
%           'MarkerFaceColor','flat', ...
%           'FaceLighting',    'gouraud',     ...
%              'AmbientStrength', 0.15);
%       title('color from Sampling');     
%     axis vis3d;
    
    
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
    C2{end+1} = sColor;
    I2{end+1} = sIndex;
%     V2{end+1} = fvs.vertices';
    
    V_all{end+1} = fv{data_idx}.Vertices';
end
V = V2';
% return;


% return;

% ground truth rotation matrices Rgt{j}*V{j} is aligned with Rgt{1}*R{1}
% Rgt = arrayfun(@(theta) angle2rotation(theta),theta,'uniformoutput',false);

% colors for each view
% clrmap = {[1 .1412 0]; [.1373 .4196 .5569]; [0 0 1]; [.8039 .6078 .1137]};
clrmap = {[1 1 0]; [1 0 1] ; [0 1 1] ; [1 0 0] ; [0 1 0] ; [0 0 1] ; [1 1 1] ; [0 0 0]};
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
marker = {'s'; 'x'; '.'; '^';'s'; 'x'; '.'; '^'};
marker = marker(1:d_size);

% initialize GMM means Xin, using random sampling of a unit sphere. Choose
% your own initialization. You may want to initialize Xin with some of the
% sets.

% set K as the 50% of the median cardinality of the views
K = ceil(0.5*median(cellfun(@(V) size(V,2),V))); 
K = 1000;

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

[R,t,X,S,a,pk,T] = jrmpc2(V,Xin,'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0);
% toc;
% tic;
% [R,t,X,S,a,pk,T] = jrmpc2_gpu(V,Xin,'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 1);
toc;

% measure and display convergency, view 1 is ommited as is the referential.
fprintf('                  ||Rgt{j} - R{j}^T*R{1}||_F                  \n');

fprintf('______________________________________________________________\n');

fprintf('Set  :'),for j=2:M,fprintf('    %d    ',j),end,fprintf('\n');

% fprintf('Error:'),for j=2:M,fprintf('  %.4f ',norm(Rgt{j}-R{j}'*R{1},'fro'));end

fprintf('\n');


% visualize the registration process, see documentation of jrmpc for T.

iter=maxNumIter;
TV_all = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V_all',T(:,1,iter),T(:,2,iter),'uniformoutput',false);

h=figure;

for iter = 1:5:maxNumIter
    % apply transformation of iteration : iter
    TV = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V,T(:,1,iter),T(:,2,iter),'uniformoutput',false);
    
    clf;
    
    hold on, grid on
    
    title(sprintf('Registration of the sets after %d iteration(s).\n',iter),'fontweight','bold','fontsize',12);
    
        
    hg2 = cellfun(@(TV,clrmap,marker,markerSize) scatter3(TV(1,:),TV(2,:),TV(3,:),markerSize.*1.5,[0 0 0],'filled'), TV, clrmap, marker, markerSize, 'UniformOutput', false);
    hg2 = cellfun(@(TV,clrmap,marker,markerSize) scatter3(TV(1,:),TV(2,:),TV(3,:),markerSize,clrmap,'filled'), TV, clrmap, marker, markerSize, 'UniformOutput', false);
    
    legend(strIdx{:});
    
%     set(gca, 'position',get(gca,'position')+[+580 0 0 0]);
    
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
axis vis3d;
saveas(h, sprintf('%s/after_reg', resultPath));
close(h);
% detect and remove "bad" centers and "unreliable" points 
[TVrefined,Xrefined,Xrem] = removePointsAndCenters(TV,X,S,a);



% figure; scatter3(Xrem(1,:)', Xrem(2,:)', Xrem(3,:)');
% figure;cellfun(@(TVrefined,clrmap,marker,mkSize) scatter3(TVrefined(1,:),TVrefined(2,:),TVrefined(3,:),mkSize,clrmap,marker),TVrefined,clrmap,marker,markerSize, 'UniformOutput', false);
    
% visualize TVrefined.
h=figure;
hold on, grid on

    title('Final registration with unreliable points removed','fontweight','bold','fontsize',12);
    
    
    hg3 = cellfun(@(TVrefined,clrmap,marker,mkSize) scatter3(TVrefined(1,:),TVrefined(2,:),TVrefined(3,:),mkSize.*1.5,[0 0 0],'filled'),TVrefined,clrmap,marker,markerSize, 'UniformOutput', false);
    hg3 = cellfun(@(TVrefined,clrmap,marker,mkSize) scatter3(TVrefined(1,:),TVrefined(2,:),TVrefined(3,:),mkSize,clrmap,'filled'),TVrefined,clrmap,marker,markerSize, 'UniformOutput', false);
    
    
    legend(strIdx{:});
    
    % use the same axes as in the registration process
%     set(gca,'XLim',XLim,'YLim',YLim,'ZLim',Zlim,'fontweight','bold','children', hg3);
    set(gca,'XLim',XLim,'YLim',YLim,'ZLim',Zlim,'fontweight','bold');
    
%     set(3,'position',get(1,'position')+[0 -510 0 0]);
    
    view([40 54]) 
hold off
saveas(h, sprintf('%s/final_remove_outlier', resultPath));
close(h);

% Visualize bad centers (orange) and good centers (blue).
h=figure;
hold on, grid on

    title('Final GMM means.','fontweight','bold','fontsize',12);
    
    scatter3(Xrefined(1,:),Xrefined(2,:),Xrefined(3,:),8,[0 .38 .67],'s');
    
    scatter3(Xrem(1,:),Xrem(2,:),Xrem(3,:),40,[1 .1412 0],'marker','x');
    
    legend('"Good" Centers','"Bad" Centers');
    
    % use the same axes as in the registration process
    set(gca,'XLim',XLim,'YLim',YLim,'ZLim',Zlim,'fontweight','bold');
    
    %set(4,'position',get(1,'position')+[+580 -510 0 0]);
    
    view([40 54])
    
hold off
saveas(h, sprintf('%s/final_gmm', resultPath));
close(h);

% return;


% sqe = @(Y,X) sum(bsxfun(@minus,permute(Y,[2 3 1]),permute(X,[3 2 1])).^2,3);
% ddd = cellfun(@(TV) sqe(TV,X),TV,'uniformoutput',false);
% ddd2 = mean([ddd{1} ; ddd{2}],1);
% ddd2 = mean([ddd{1}],1);
% 
% figure;scatter3(X(1,:),X(2,:),X(3,:),40,ddd2,'filled'); colorbar;
% figure;scatter3(X(1,:),X(2,:),X(3,:),40,pk,'filled'); colorbar;
% figure;scatter3(X(1,:),X(2,:),X(3,:),40,S,'filled'); colorbar;


A = zeros(K, d_size);
%A2 = zeros(K, d_size);
%A3 = zeros(K, d_size);
%A4 = zeros(K, d_size);%zeros(K, 1);

maxIdx = {};

for ii=1:d_size
    [m, i] = max(a{ii}');
    maxIdx{ii} = i;
    gmmCount = accumarray(i', 1);
    A(1:length(gmmCount), ii) = gmmCount;

    %A = [A accumarray(i', 1)];
    %A2 = [A2 sum(a{ii})'];
    %A3 = [A3 mean(a{ii})'];

    [sortedValues,sortIndex] = sort(a{ii}','descend'); 


    firstIdx = sortIndex(1, :)';
    firstVal = sortedValues(1, :)';
    secondIdx = sortIndex(2, :)';
    secondVal = sortedValues(2, :)';

    % KKK = accumarray(firstIdx, firstVal);
    % KKK2 = accumarray(secondIdx, secondVal);
    % 
    % A4 = [A4 KKK+KKK2];
end
% AA = bsxfun(@rdivide, A, max(A')');
% figure; plot(AA');
h=figure; plot(A'); title(sprintf('Group size of each GMM (#K : %d)', K));
saveas(h, sprintf('%s/groupassign_k%d', resultPath, K));
close(h);
% figure; plot(A2');
% figure; plot(A3');
% figure; plot(A4');

% gmm_score = AA;
% 
% cluster = [];
% for i=1:K
%     cluster = [cluster vl_kmeans(gmm_score(i,:), 2)];
% end
% figure; plot(cluster)


% return;


for id=1:d_size
    fv3.Faces = fv{id}.Faces;
    fv3.Vertices = fv{id}.Vertices;

    sColor = C2{id};

    mIdx = maxIdx{id};
    sRMean2 = accumarray(mIdx', sColor(1,:)', [], @mean);
    sGMean2 = accumarray(mIdx', sColor(2,:)', [], @mean);
    sBMean2 = accumarray(mIdx', sColor(3,:)', [], @mean);

    GMM_color = [sRMean2 sGMean2 sBMean2];
    sampleColor = GMM_color(mIdx, :)';

    sIndex = I2{id};
    fv3.FaceVertexCData = sampleColor(:,sIndex)';

    h = figure;
    patch(fv3,'FaceColor','flat','EdgeColor','flat',...
          'MarkerFaceColor','flat', ...
          'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);
    title('color from GMM');     
    axis vis3d;   
    saveas(h, sprintf('%s/color_GMM_%d', resultPath, id));
    close(h);

    
    
    rGMM_color = bsxfun(@plus, rand(K, 3)./2, [0.3 0.3 0.3]);
    sampleColor = rGMM_color(mIdx, :)';
    sIndex = I2{id};
    fv3.FaceVertexCData = sampleColor(:,sIndex)';

    h = figure;
    patch(fv3,'FaceColor','flat','EdgeColor','flat',...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);
    title('color from GMM');     
    camlight('headlight');
    material('dull');
    axis vis3d;   
    saveas(h, sprintf('%s/rcolor_GMM_%d', resultPath, id));
    close(h);
    

    fv2.Faces =  fv{id}.Faces;
    fv2.Vertices = fv{id}.Vertices;
    sColor = C2{id};
    fv2.FaceVertexCData = sColor(:,I2{id})';
    % 
    h = figure;
    patch(fv2,'FaceColor','flat','EdgeColor','flat',...
          'MarkerFaceColor','flat', ...
          'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);
      title('color from Sampling');     
    axis vis3d;
    saveas(h, sprintf('%s/color_sampling_%d.fig', resultPath, id));
    close(h);
end

AX = [A' ;X.*sqrt(size(A,2)).*mean(mean(A)).*1];
labels = '_AX2';
% [centers, assignments] = vl_kmeans(A', 16);
[centers, assignments] = vl_kmeans(AX, 16);

A2 = A';

rmax = 4;
cmax = 4;

h=figure;
pIdx = 1;
for rr=1:rmax
    for cc=1:cmax
        subplot(rmax, cmax, pIdx);
        hold on;
        A2sub = A2(:, find(assignments == pIdx));
        [m, i] = max(mean(A2sub'));
        plot(A2(:, find(assignments == pIdx)));
        plot(mean(A2sub'), 'k', 'lineWidth', 2);
        title(sprintf('%d-th, size : %d, var:%.2f', pIdx, sum(assignments == pIdx), mean(var(A2sub'))));
        hold off;
        pIdx = pIdx+1;
    end
end
saveas(h, sprintf('%s/groupassign_kmeans%s', resultPath, labels));
close(h);

h=figure;
pIdx = 1;
for rr=1:rmax
    for cc=1:cmax
        subplot(rmax, cmax, pIdx);
        hold on;
        A2sub = A2(:, find(assignments == pIdx));
        [m, i] = max(mean(A2sub'));
        
        norm_a2 = bsxfun(@minus, A2(:, find(assignments == pIdx)), mean(A2sub));
        norm_a2 = norm_a2 ./ sqrt(var(mean(A2sub)));
%         norm_a2 = bsxfun(@rdivide, norm_a2, sqrt(var(A2sub'))');
        plot(norm_a2);
        plot(mean(norm_a2'), 'g', 'lineWidth', 2);
%         plot(mean(A2sub'), 'k', 'lineWidth', 2);
        norm_var = mean(var(norm_a2'));
        
        plot(ones(size(norm_a2, 1), 1).*norm_var.*0.5, 'k', 'lineWidth', 2);
        
        title(sprintf('%d-th, size : %d, norm-var:%.2f', pIdx, sum(assignments == pIdx), norm_var));
        hold off;
        pIdx = pIdx+1;
    end
end
saveas(h, sprintf('%s/groupassign_kmeans_norm%s', resultPath, labels));
close(h);



h=figure;
pIdx = 1;
for rr=1:rmax
    for cc=1:cmax
        subplot(rmax, cmax, pIdx);
        Xg = X(:,find(assignments == pIdx));
        Xother = X(:,find(assignments ~= pIdx));
        hold on;
        scatter3(Xg(1,:),Xg(2,:),Xg(3,:),8,[0 .38 .67],'s');
        scatter3(Xother(1,:),Xother(2,:),Xother(3,:),8,[.67 .1 0],'s');
        hold off;
%         plot(A2(:, find(assignments == pIdx)));
        title(sprintf('size : %d', sum(assignments == pIdx)));
        pIdx = pIdx+1;
    end
end
saveas(h, sprintf('%s/groupassign_kmeans_gmm%s', resultPath, labels));
close(h);


id=1;
%GMM group index
gmmk = 1;
%included GMM index in goup GMMK

% plotGMMK( id, gmmk, assignments, maxIdx, fv, I2 );
id=8;
gmmk=14;
h=plotGMMK( id, gmmk, assignments, maxIdx, fv, I2 );
title(sprintf('%d-th set, gmm group=%d', id, gmmk));
saveas(h, sprintf('%s/marked_set_%d_gmmk_%d%s', resultPath, id, gmmk, labels));
close(h);

id=3;
gmmk=13;
h=plotGMMK( id, gmmk, assignments, maxIdx, fv, I2 );
title(sprintf('%d-th set, gmm group=%d', id, gmmk));
saveas(h, sprintf('%s/marked_set_%d_gmmk_%d%s', resultPath, id, gmmk, labels));
close(h);


id=1;
gmmk=15;
h=plotGMMK( id, gmmk, assignments, maxIdx, fv, I2 );
title(sprintf('%d-th set, gmm group=%d', id, gmmk));
saveas(h, sprintf('%s/marked_set_%d_gmmk_%d%s', resultPath, id, gmmk, labels));
close(h);

id=6;
gmmk=10;
h=plotGMMK( id, gmmk, assignments, maxIdx, fv, I2 );
title(sprintf('%d-th set, gmm group=%d', id, gmmk));
saveas(h, sprintf('%s/marked_set_%d_gmmk_%d%s', resultPath, id, gmmk, labels));
close(h);

id=5;
gmmk=15;
h=plotGMMK( id, gmmk, assignments, maxIdx, fv, I2, 2 );
title(sprintf('%d-th set, gmm group=%d', id, gmmk));
saveas(h, sprintf('%s/marked_set_%d_gmmk_%d%s', resultPath, id, gmmk, labels));
close(h);

% gmmk_gmmidx = find(assignments==gmmk);
% %included sample index
% gmmk_sampleIdx = find(ismember(maxIdx{id}, gmmk_gmmidx));
% %included original Point cloud index
% gmmk_Idx = find(ismember(I2{id}, gmmk_sampleIdx));
% 
% fv2.Faces =  fv{id}.Faces;
% fv2.Vertices = fv{id}.Vertices;
% sColor = C2{id};
% fv2.FaceVertexCData = fv{id}.FaceVertexCData;%sColor(:,I2{id})';
% fv2.FaceVertexCData(gmmk_Idx,:) = bsxfun(@plus,fv2.FaceVertexCData(gmmk_Idx,:),[0.5 0 0]);
% % 
% h = figure;
% patch(fv2,'FaceColor','flat','EdgeColor','flat',...
%       'MarkerFaceColor','flat', ...
%       'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);
%   title('color from Sampling');     
% axis vis3d;
