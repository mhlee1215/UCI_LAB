clc
close all
clear all

addpath(genpath('libs'));
isDebug = true;


[v,f,c] = loadData3D('~/Desktop/3dScene/1Bedroom/77_labels.obj');
[v2,f2,c2] = loadData3D('~/Desktop/3dScene/small/chair/chair_1.ply');
[v3,f3,c3] = loadData3D('~/Desktop/3dScene/small/chair/table3_1.ply');

chairMat = read4x4Mats('syn_data/chair_locs.txt');
tableMat = read4x4Mats('syn_data/table_locs.txt');

f2 = f2 + size(v, 1);
f3 = f3 + size(v, 1) + size(v2, 1);
f_all = [f ; f2 ; f3];



for pi=1%:length(chairMat)

    v2_t = bsxfun(@plus,chairMat{pi}.R*v2',chairMat{pi}.t)';
    v3_t = bsxfun(@plus,tableMat{pi}.R*v3',tableMat{pi}.t)';


    v_all = [v ; v2_t ; v3_t];
    resultFilePath = sprintf('bedroom_syn_parametric_test_%d.klg', pi);

    % figure;
    % fv2.Faces = f_all;
    % fv2.Vertices = v_all;
    % fv2.FaceVertexCData = zeros(size(fv2.Vertices));%sColor(:,I2{id})';
    % fv2.FaceVertexCData = bsxfun(@plus, fv2.FaceVertexCData, [0.7 0.7 0.9]);
    % % fv2.FaceVertexCData(gmmk_Idx,:) = bsxfun(@plus,bsxfun(@times, fv2.FaceVertexCData(gmmk_Idx,:), [0 0 0]),[0.9 0.4 0.4]);
    % patch(fv2, ...
    %      'FaceColor','flat','EdgeColor','flat', ...
    %      'EdgeColor',       'none',        ...
    %      'FaceLighting',    'gouraud',     ...
    %      'AmbientStrength', 0.15);
    % camlight('headlight');
    % material('dull');
    % axis equal;



  
    % dataPath = '~/3d_dataset/1Bedroom/77_labels_2.obj';
    poseFilePath = 'syn_path/out.txt';
    ScreenSizeV = [480;640];



    fprintf(' Data loading...\n');
    % obj = loadawobj(dataPath);
    % vertices = obj.v';
    % faces = obj.f3';

    vertices = v_all;
    faces = f_all;


    %Load synthetic pose set

    poseSet = EF_load_path(poseFilePath, vertices);

    delete(resultFilePath);
    logFile = fopen(resultFilePath,'wb+');
    maxFrame = length(poseSet);
    fwrite(logFile, maxFrame, 'int32');
    if isDebug
        close all; figure;     
    end
    for ii=1:maxFrame
        sprintf('%d/%d ...\n', ii, maxFrame)
        TT = poseSet{ii}.t;
        RR = poseSet{ii}.R;

        cam.fcV = [500 ; 500];
        cam.TcV = (TT).* [-1 1 -1]';%RMat*(T+TMat);%RMat*TMat+T;%RMat*(T+TMat);%RMat*T+TMat .* [1 -1 -1]';
        cam.RcM = RR .* [-1 -1 -1; 1 1 1; -1 -1 -1];%RMat*R;%*RMat;%RMat*R .* [1 1 1; -1 -1 -1; -1 -1 -1];
        cam.ccV = [ScreenSizeV(2)./2 ScreenSizeV(1)./2]';

        invertedDepth = false;

        zoomFactor = [1 1];
        ZNearFarV = [0.1 10]';
        [D, CameraCoordT] = RenderDepthMesh(faces, vertices, cam, ScreenSizeV, ZNearFarV, zoomFactor, invertedDepth);
        if isDebug
            imagesc(CameraCoordT(:,:,3)); colorbar;
            pause(0.01);
        end

        DD = CameraCoordT(:,:,3)./10;
        DepthImageM = int16(DD.*4096);
        fwrite(logFile, 1000, 'int64');
        fwrite(logFile, prod(ScreenSizeV)*2, 'int32');
        fwrite(logFile, 0, 'int32');
        fwrite(logFile, repmat(DepthImageM', 1, 1), 'ushort');
    end
    fclose(logFile);

end
