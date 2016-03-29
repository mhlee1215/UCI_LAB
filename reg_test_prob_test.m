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
data_set = 1:10;%[1 2 3 4];
d_size = length(data_set);

resultPath = 'results/com_with_without_vis';
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
    
    dataRoot = 'E:\lab\work\2015_4d_seg\codes\data\desk_color3';
    [tri, pts, data, comments] = ply_read(sprintf('%s/MeshedReconstruction_%d.ply', dataRoot, data_idx), 'tri');
    
    fprintf('end\n');
    
    
    
    
    v = [data.vertex.x data.vertex.y data.vertex.z];
    f = tri';
    c = [data.vertex.red data.vertex.green data.vertex.blue];

    fv{data_idx}.Faces = f;
    fv{data_idx}.Vertices = v;
    fv{data_idx}.FaceVertexCData = c./255;

    [sData, sIndex, sColor] = uniformSubSample(fv{data_idx}.Vertices, 70, fv{data_idx}.FaceVertexCData);
    
%     if data_idx == 1
%         normVar = sqrt(var(sData')');
%     end
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
clrmap = {[1 1 0]; [1 0 1] ; [0 1 1] ; [1 0 0] ; [0 1 0] ; [0 0 1] ; [1 1 1] ; [0 0 0] ; [0.3 0.3 0.3] ; [.5 0.9 0.9] ; [0.8 1 1]; [0 0 0]};
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
marker = {'s'; 'x'; '.'; '^';'s'; 'x'; '.'; '^';'s'; 'x'; '.'; '^'};
marker = marker(1:d_size);

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

[R,t,X,S,a,pk,T,TAssigned, TXQ] = jrmpc_soft(V,Xin,'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0);
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

%  figure;
%  for i=length(T)
%      curX = TXQ{1,i};
%      curQ = TXQ{2,i};
%      for j=1:size(curX, 2)
%         h1 = plot_gaussian_ellipsoid(curX(:,j), [curQ(j) 0 0 ; 0 curQ(j) 0 ; 0 0 curQ(j)]./max(curQ), curQ(j)/max(curQ)*0.03);
%      end
%  end
%  
%  set(h1,'facealpha',0.6);
%  view(129,36); set(gca,'proj','perspective'); grid on; 
%  grid on; axis equal; axis tight;

params.type = 2;
params.Assigned = TAssigned;
params.K = K;
params.interval = 3;
params.marker = marker;
params.markerSize = markerSize;
params.clrmap = clrmap;
params.strIdx = strIdx;

h = figure;
[TV] = drawTransformation(V, T, params);
axis equal;

% return;
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
paramsPlotGMMK.maxIdx = TAssigned2(:,end);
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
paramsPlotGMMK.maxIdx = TAssigned5(:,end);
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

