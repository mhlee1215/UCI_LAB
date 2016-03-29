clc
close all
clear all
% g = gpuDevice(1);
% reset(g);

addpath(genpath('libs'));
addpath(genpath('libs2'));
run('libs2/vlfeat-0.9.20-bin/vlfeat-0.9.20/toolbox/vl_setup');
data_set = [1:6 8:14];%[1 2 3 4];
d_size = length(data_set);

resultPath = 'results/syn_data3';
mkdir(resultPath);



fprintf(' Data loading...\n');

obj = loadawobj('1_labels.obj');
vertices = obj.v';
faces = obj.f3';



ScreenSizeV = [480 ; 640];
cam.fcV = [525 ; 525];
cam.TcV = [0 2 3]';%.*[1 -1 -1]';%(TT).* [-1 1 -1]';%RMat*(T+TMat);%RMat*TMat+T;%RMat*(T+TMat);%RMat*T+TMat .* [1 -1 -1]';
% theta = (35*pi)/180;
% cam.RcM = [1 0 0 ; 0 cos(theta) -sin(theta) ; 0 sin(theta) cos(theta)].*[1 1 1; -1 -1 -1; -1 -1 -1];
cam.RcM = [1 0 0 ; 0 1 0 ; 0 0 1].*[1 1 1; -1 -1 -1; -1 -1 -1];


% cam.RcM = [1 0 0 ; 0 1 0 ; 0 0 1].*[1 1 1; -1 -1 -1; -1 -1 -1];%.*[-1 -1 -1; 1 1 1; -1 -1 -1];% .* [-1 0 0 ; 0 1 0 ; 0 0 -1];%.* [1 1 1; -1 -1 -1; -1 -1 -1];%RR .* [-1 -1 -1; 1 1 1; -1 -1 -1];%RMat*R;%*RMat;%RMat*R .* [1 1 1; -1 -1 -1; -1 -1 -1];
cam.ccV = [ScreenSizeV(2)./2 ScreenSizeV(1)./2]';

%     vertices = fv{j}.small.vertices;
%     vertices = bsxfun(@plus, vertices', TMat)';
%     vertices = (RMat*vertices')';%bsxfun(@plus, RMat*vertices', TMat)';
%     faces = fv{j}.small.faces;
invertedDepth = false;

zoomFactor = [1 ; 1];
ZNearFarV = [0.1 10]';%[0; 0];
[DepthImageM, CameraCoordT] = RenderDepthMesh(faces, vertices, cam, ScreenSizeV, ZNearFarV, zoomFactor, invertedDepth);
figure; imagesc(CameraCoordT(:,:,3)); colorbar;
% 
% return;









% h=figure; hold on;
% fv_alligned.Faces = obj.f3';%delaunay(vertices(:,1),vertices(:,2));%fv_small.faces;
% fv_alligned.Vertices = obj.v';%fv_small.vertices;
% patch(fv_alligned,'FaceColor',       [0.8 0.8 1.0], ...
%          'EdgeColor',       'none',        ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);
% camlight('headlight');
% material('dull');
% axis equal;

%----------------------------------------------------
% [tri, pts, data, comments] = ply_read('livingroom.ply', 'tri');
% vertices = [data.vertex.x data.vertex.y data.vertex.z];
% faces = tri';

poseFilePath = 'out_2.txt';
poseFileID = fopen(poseFilePath, 'r');
poseMat = fscanf(poseFileID, '%f');
poseSet = {};
% imW = 3;
% imH = 3;
for pii = 1:size(poseMat, 1)/19
    idx = (pii-1)*19+1;
    pM = reshape(poseMat(idx+3:idx+18)', 4, 4)';
    cameraMat.R = pM(1:3, 1:3);
    cameraMat.t = pM(1:3, end);
    cameraMat.imSize = [0 0];%[imH imW];
    cameraMat.f = 1;
    poseSet{end+1} = cameraMat;
end

% [vIdx, inVIdx, inVvec] = HPR(vertices,[3 3 3],1.2);
% [c_iv, c_if] = cropFV(vertices', faces, inVIdx);
% vertices2 = c_iv';
% faces2 = c_if;

% figure; scatter3(vertices(:,1), vertices(:,2), vertices(:,3), 0.5, 'filled');

% faces = delaunay(vertices(:,1),vertices(:,2));
%write_ply(vertices, faces, 'delunay.ply');

%noOfpoints = length(mesh.vertices);
%noOfpolygons = length(mesh.triangles);
% mesh.vertices = vertices;
% mesh.triangles = faces;
% makePly(mesh, 'delunay2.ply');

% h=figure; hold on;
% fv_alligned.Faces = faces;%delaunay(vertices(:,1),vertices(:,2));%fv_small.faces;
% fv_alligned.Vertices = vertices;%fv_small.vertices;
% patch(fv_alligned,'FaceColor',       [0.8 0.8 1.0], ...
%          'EdgeColor',       'none',        ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);
% camlight('headlight');
% material('dull');
% axis equal;
% hold on;
% return;
% for i=1:length(poseSet)
%     cameraMat = poseSet{i};
%     camPos = -cameraMat.R'*cameraMat.t;
%     scatter3(camPos(1), camPos(2), camPos(3), 0.3);
% end
% for i=1:100:length(poseSet)
%     cameraMat = poseSet{i};
% %    pM(:,end)
% %     scatter3(cameraMat.t(1), cameraMat.t(2), cameraMat.t(3), 0.1);
% %     P = projectImg2World(pM, [0 0]', 0.1);
%     
%     
%     cameraMat.R
%     Cp = -cameraMat.R'*cameraMat.t;
%     depth = 0.5;
%     C1 = projectImg2World(cameraMat, [-imH -imW]', depth);
%     C2 = projectImg2World(cameraMat, [-imH imW]', depth);
%     C3 = projectImg2World(cameraMat, [imH imW]', depth);
%     C4 = projectImg2World(cameraMat, [imH -imW]', depth);
%     PScene = [C1 ; C2 ; C3 ; C4 ; C1];
% 
%     for i=1:size(PScene, 1)-1
%         line([Cp(1) PScene(i,1)]', [Cp(2) PScene(i,2)]', [Cp(3) PScene(i,3)]', 'linewidth', 1);    
%         line([PScene(i,1) PScene(i+1,1)]', [PScene(i,2) PScene(i+1,2)]', [PScene(i,3) PScene(i+1,3)]', 'linewidth', 1);    
%     end
%     
% end
% axis equal;
% saveas(h, sprintf('%s/trajectory_1.fig', resultPath));
% close(h);

% [vertices_dense, faces_dense] = linearSubdivision(vertices', faces');


ScreenSizeV = round([atan(62.7)*525 ; atan(62.7)*525]);

% cam.TcV = [0 2 1]';%.*[1 -1 -1]';%(TT).* [-1 1 -1]';%RMat*(T+TMat);%RMat*TMat+T;%RMat*(T+TMat);%RMat*T+TMat .* [1 -1 -1]';
% theta = (35*pi)/180;
% % cam.RcM = [1 0 0 ; 0 cos(theta) -sin(theta) ; 0 sin(theta) cos(theta)].*[1 1 1; -1 -1 -1; -1 -1 -1];
% cam.RcM = [1 0 0 ; 0 1 0 ; 0 0 1].*[1 1 1; -1 -1 -1; -1 -1 -1];
    
    
delete('livingroom_syn.klg');
logFile = fopen('livingroom_syn.klg','wb+');

maxFrame = length(poseSet);
fwrite(logFile, maxFrame, 'int32');
%   figure;     
for ii=1:20%length(poseSet)
    ii
    TT = poseSet{ii}.t;
    RR = poseSet{ii}.R;

    cam.fcV = [525 ; 525];
    cam.TcV = (TT).* [-1 1 -1]';%RMat*(T+TMat);%RMat*TMat+T;%RMat*(T+TMat);%RMat*T+TMat .* [1 -1 -1]';
    cam.RcM = RR .* [-1 -1 -1; 1 1 1; -1 -1 -1];%RMat*R;%*RMat;%RMat*R .* [1 1 1; -1 -1 -1; -1 -1 -1];
    cam.ccV = ScreenSizeV ./2;

    %     vertices = fv{j}.small.vertices;
    %     vertices = bsxfun(@plus, vertices', TMat)';
    %     vertices = (RMat*vertices')';%bsxfun(@plus, RMat*vertices', TMat)';
    %     faces = fv{j}.small.faces;
    invertedDepth = false;

    zoomFactor = [1 1];
    ZNearFarV = [0.1 10]';%[0; 0];
    [D, CameraCoordT] = RenderDepthMesh(faces, vertices, cam, ScreenSizeV, ZNearFarV, zoomFactor, invertedDepth);
    %imagesc(CameraCoordT(:,:,3)); colorbar;
   
    DepthImageM = CameraCoordT(:,:,3);
    
   DepthImageM = int16(DepthImageM.*4096);
   
    fwrite(logFile, 1000, 'int64');
    fwrite(logFile, prod(ScreenSizeV)*2, 'int32');
    fwrite(logFile, 0, 'int32');
    fwrite(logFile, repmat(DepthImageM, 1, 1), 'int16');
%     fwrite(logFile, encodedImage->data.ptr, 'float32');
    
     %pause(0.01);
end
    
fclose(logFile);



% 
% testF = fopen('/home/mhlee/Downloads/ElasticFusion-master/ElasticFusion-master/GUI/build/2016-02-03.00.klg', 'rb+');
% 
% fread(testF, 1, 'int64')
% fread(testF, 1, 'int32')
% fread(testF, 1, 'int32')
% 
























return;









for ii=1:100:length(poseSet)
    cameraMat = poseSet{ii};
    Cp = -cameraMat.R'*cameraMat.t;
%    pM(:,end)
%     scatter3(cameraMat.t(1), cameraMat.t(2), cameraMat.t(3), 0.1);
%     P = projectImg2World(pM, [0 0]', 0.1);
    
    v2 = -cameraMat.R'*(bsxfun(@plus, vertices', cameraMat.t));
    vt = v2;
    
    positive = (v2(3,:) < 0);
    v2(1,:) = v2(1,:) ./ v2(3,:);
    v2(2,:) = v2(2,:) ./ v2(3,:);
%     v2(3,:) = v2(3,:) ./ v2(3,:);
    
    %check whetehr 3D point is visible or not
    %[vIdx, inVIdx, inVvec] = HPR(vertices,Cp',1.8);
    %check whether point is correctly projected into 2D screen
    valid_idx = find(positive .* (v2(1,:)> -imH).*(v2(1,:) < imH) .* (v2(2,:) > -imW) .* (v2(2,:) < imW));
%     valid_idx = find((1-inVvec)'.*positive .* (v2(1,:)> -imH).*(v2(1,:) < imH) .* (v2(2,:) > -imW) .* (v2(2,:) < imW));
    [c_v, c_f] = cropFV(vertices', faces, valid_idx);
%     valid_face_idx = find(ismember(faces(:,1), valid_idx).*ismember(faces(:,2), valid_idx).*ismember(faces(:,3), valid_idx));
%     v3 = vertices';
%     v_crop = v3(:,valid_idx);
% k = valid_idx;
%     v = 1:length(valid_idx);
%     map = containers.Map(k, v);
% 
%     f1_cvted = arrayfun(@(x) map(x), faces(valid_face_idx, 1));
%     f2_cvted = arrayfun(@(x) map(x), faces(valid_face_idx, 2));
%     f3_cvted = arrayfun(@(x) map(x), faces(valid_face_idx, 3));

    fv_crop.vertices = c_v;
    fv_crop.faces = c_f;
    h=figure; hold on;
    scatter3(c_v(1,:), c_v(2,:), c_v(3,:), 0.1);
    
%     cameraMat.R
    
    depth = 0.5;
    C1 = projectImg2World(cameraMat, [-imH -imW]', depth);
    C2 = projectImg2World(cameraMat, [-imH imW]', depth);
    C3 = projectImg2World(cameraMat, [imH imW]', depth);
    C4 = projectImg2World(cameraMat, [imH -imW]', depth);
    PScene = [C1 ; C2 ; C3 ; C4 ; C1];

    for i=1:size(PScene, 1)-1
        line([Cp(1) PScene(i,1)]', [Cp(2) PScene(i,2)]', [Cp(3) PScene(i,3)]', 'linewidth', 1);    
        line([PScene(i,1) PScene(i+1,1)]', [PScene(i,2) PScene(i+1,2)]', [PScene(i,3) PScene(i+1,3)]', 'linewidth', 1);    
    end
    
    axis equal;
    
    saveas(h, sprintf('%s/view_%d.fig', resultPath, ii));
    close(h);

%     h=figure; fv_alligned.Faces = faces;%fv_crop.faces;%fv_small.faces;
%     fv_alligned.Vertices = vt';%fv_crop.vertices';%fv_small.vertices;
%     patch(fv_alligned,'FaceColor',       [0.8 0.8 1.0], ...
%              'EdgeColor',       'none',        ...
%              'FaceLighting',    'gouraud',     ...
%              'AmbientStrength', 0.15);
%     camlight('headlight');
%     material('dull');
%     axis equal;
%     hold on;
    
    [c_v, c_f] = cropFV(vt, faces, valid_idx);
    h=figure; fv_alligned.Faces = c_f;%fv_crop.faces;%fv_small.faces;
    fv_alligned.Vertices = c_v';%fv_crop.vertices';%fv_small.vertices;
    patch(fv_alligned,'FaceColor',       [0.8 0.8 1.0], ...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);
    camlight('headlight');
    material('dull');
    axis equal;
    hold on;
    saveas(h, sprintf('%s/model_%d.fig', resultPath, ii));
    close(h);
end


return;


% fprintf('Generate small version..');
% fv_small.vertices = v;
% fv_small.faces = f;
% fv_small = reducepatch(fv_small, 0.01);
% fprintf('End\n');
% fv{data_idx}.small = fv_small;

 figure; fv_alligned.Faces = faces2;%fv_small.faces;
    fv_alligned.Vertices = v;%fv_small.vertices;
    patch(fv_alligned,'FaceColor',       [0.8 0.8 1.0], ...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);
    camlight('headlight');
    material('dull');
    axis equal;
    hold on;
    
   


return;


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
marker = {'s'; 'x'; '.'; '^';'s'; 'x'; '.'; '^';'s'; 'x'; '.'; '^';'s'...
    ; 'x'; '.'; '^';'s'; 'x'; '.'; '^';'s'; 'x'; '.'; '^'};
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

TV_small_all = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V_all',T(:,1,iter),T(:,2,iter),'uniformoutput',false);

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
params.view = [];%[40 54];
params.pause = .0;
h = figure;
[TV] = drawTransformation(V, T, params);
axis equal;

% return;



% for j=1%:M    
%     RMat = T{j,1,end};
%     TMat = T{j,2,end};
%     
% %     poseMat = cell2mat(Pose_Set{j})';
% 
%     ScreenSizeV = [200 200]';
% 
%     
% %     C2 = [[eye(3) TMat];0 0 0 1];
%     C2 = [[eye(3) zeros(3,1)];0 0 0 1];
%     
%     Ccell = PoseMat_Set{j};
%     C = Ccell{1};
%     C = C*C2;
%     C = C(1:3, :);
%     M = C(:, 1:3);
%     M = M;
% %     [K, RR] = rq(M);
%     [K, RR] = rq2(M);
% %     RR = inv(K)*RR;
% %     K=K*K;
%     TT = M\(C(:,end));
% 
% 
%     cam.fcV = [K(1,1) K(2,2)]';
%     cam.TcV = (TT).* [-1 1 -1]';%RMat*(T+TMat);%RMat*TMat+T;%RMat*(T+TMat);%RMat*T+TMat .* [1 -1 -1]';
%     cam.RcM = RR .* [-1 -1 -1; 1 1 1; -1 -1 -1];%RMat*R;%*RMat;%RMat*R .* [1 1 1; -1 -1 -1; -1 -1 -1];
%     cam.ccV = ScreenSizeV ./2;
%     
%     vertices = fv{j}.small.vertices;
% %     vertices = bsxfun(@plus, vertices', TMat)';
% %     vertices = (RMat*vertices')';%bsxfun(@plus, RMat*vertices', TMat)';
%     faces = fv{j}.small.faces;
%     invertedDepth = false;
%     
%     zoomFactor = 1;
%     ZNearFarV = [0 200]';
%     [DepthImageM, CameraCoordT] = RenderDepthMesh(faces, vertices, cam, ScreenSizeV, ZNearFarV, zoomFactor, invertedDepth);
%     D = CameraCoordT(:,:,3);
% %     figure; imagesc(CameraCoordT(:,:,3)); colorbar;
%     figure; imagesc(DepthImageM); colorbar;
%     
%     figure; fv_alligned.Faces = faces;
%     fv_alligned.Vertices = vertices;
%     patch(fv_alligned,'FaceColor',       [0.8 0.8 1.0], ...
%              'EdgeColor',       'none',        ...
%              'FaceLighting',    'gouraud',     ...
%              'AmbientStrength', 0.15);
%     camlight('headlight');
%     material('dull');
%     axis equal;
%     hold on;
%     scatter3(TT(1)+TMat(1), TT(2)++TMat(2), TT(3)++TMat(3));
%     scatter3(cam.TcV(1), cam.TcV(2), cam.TcV(3), 'k', 'filled');
% 
%     poseMat = bsxfun(@plus, RMat*poseMat', TMat)';
% end


% return;

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

