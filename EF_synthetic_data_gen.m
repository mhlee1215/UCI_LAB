% clc
% close all
% clear all

addpath(genpath('libs'));
isDebug = false;

colorType = 'color'; %'color'
randomColoring = true;
dt = 0.003;
p = [ -2.0 1.5 0;
%         -2 1.5 3;
        -0.0 1.5 2;
        2.0 1.5 -0 ;
        -0.0 1.5 -2
%         -2.5 1.5 -1;
%         -0.5 1.5 2;
%         -0.5 1.5 0;
        ];
p = [p ; p(1:2,:)];

params = {};

params{end+1}.category = 'Bedroom';
params{end}.name = 'bedroom_wenfagx';
params{end}.offset = [0.5 0 0];
params{end}.scale = [0.2 1 0.2];

params{end+1}.category = 'Bedroom';
params{end}.name = '77_labels';
params{end}.offset = [0 0 0];
params{end}.scale = [0.2 1 0.2];


[v2,f2,c2] = loadData3D('~/Desktop/3dScene/small/chair/chair_1.ply');
c2 = ones(size(c2));
c2 = bsxfun(@times, c2, [1 0 0]);
[v3,f3,c3] = loadData3D('~/Desktop/3dScene/small/chair/table3_1.ply');
c3 = ones(size(c3));
c3 = bsxfun(@times, c3, [0 1 0]);

% drawMesh(v2, f2, c2);

for dIdx=2%1:length(params)

    % [v,f,c] = loadData3D('~/Desktop/3dScene/1Bedroom/77_labels.obj');

    category = params{dIdx}.category;%'Bedroom';
    % name = '77_labels';
    name = params{dIdx}.name;%'bedroom_wenfagx';
    [v,f,c, g] = loadData3D(sprintf('~/Desktop/3dScene/1%s/%s.obj', category, name));
    if isempty(c)
        c = rand(size(v));
    end
    
    %random coloring
    if randomColoring
%         [clustCent,point2cluster,clustMembsCell] = ...
%            MeanShiftCluster([v]', 0.5);

        for ii=1:length(g)-3   
            c(g(ii):g(ii+1), 1) = rand()*0.8;
            c(g(ii):g(ii+1), 2) = rand()*0.8;
            c(g(ii):g(ii+1), 3) = rand()*0.8;
        end
    end
   

    chairMat = read4x4Mats('syn_data/chair_locs.txt');
    tableMat = read4x4Mats('syn_data/table_locs.txt');

    f22 = f2 + size(v, 1);
    f33 = f3 + size(v, 1) + size(v2, 1);
    f_all = [f ; f22 ; f33];
    c_all = [c ; c2 ; c3];


    for type = 2
    for pi=1:2%length(chairMat)

        v2_t = bsxfun(@plus,chairMat{pi}.R*v2',chairMat{pi}.t)';
        v3_t = bsxfun(@plus,tableMat{pi}.R*v3',tableMat{pi}.t)';


        v_all = [v ; v2_t ; v3_t];
        resultFilePath = sprintf('./data/syn/syn_parametric_%s_%s_%s_%.4f_%d_%d.klg', colorType, category, name, dt, type, pi);

        ScreenSizeV = [480;640];



        fprintf(' Data loading...\n');
      
        vertices = v_all;
        faces = f_all;
        color = (c_all);
        
%         drawMesh(vertices, faces, double(color));

        %Load synthetic pose set
    %     poseSet = EF_load_path(poseFilePath, vertices);
    % 
    
        offset = params{dIdx}.offset;%[0 0 0]; %for beadroom.
    %     offset = [1.5 3 0];
        scale = params{dIdx}.scale;%[0.5 0.5 0.5];
        p2 = bsxfun(@plus, p, offset);
        p2 = bsxfun(@times, p2, scale);


        %close loop
    %     p = [p ; p(1,:)];
        pc = [0 1.6 0];
        if isDebug
            poseSet = syn_path_gen(p2, dt, vertices(1:1:end, :), type, pc);
%             view(-146, 5);
            view(-180, 0);
            view(-160, -14);
        else
            poseSet = syn_path_gen(p2, dt, [], type, pc);
        end
    %     view(151, -7)

    

        
        maxFrame = length(poseSet);
        if isDebug
    %         close all; 
            h=figure;     
        else
            delete(resultFilePath);
            logFile = fopen(resultFilePath,'wb+');   
            fwrite(logFile, maxFrame, 'int32');    
        end
        for ii=1:maxFrame
            sprintf('%d/%d ...\n', ii, maxFrame)
            TT = poseSet{ii}.t;
            RR = poseSet{ii}.R;

            cam.fcV = [500 ; 500];
            cam.TcV = (TT).* [-1 1 1]';%RMat*(T+TMat);%RMat*TMat+T;%RMat*(T+TMat);%RMat*T+TMat .* [1 -1 -1]';
            cam.RcM = RR.* [-1 -1 -1; 1 1 1; 1 1 1];%.* [-1 -1 -1; 1 1 1; -1 -1 -1];%RMat*R;%*RMat;%RMat*R .* [1 1 1; -1 -1 -1; -1 -1 -1];
            cam.ccV = [ScreenSizeV(2)./2 ScreenSizeV(1)./2]';

            invertedDepth = false;

            zoomFactor = [1 1];
            ZNearFarV = [0.1 10]';
            [D, CameraCoordT] = RenderDepthMesh(faces, vertices, cam, ScreenSizeV, ZNearFarV, zoomFactor, invertedDepth);
            
            DD = CameraCoordT(:,:,3)./10;
            DepthImageM = int16(DD.*4096);

            if type == 2
                DepthImageM = flipud(DepthImageM);
            end
            
            if strcmp(colorType, 'color') == 1
                ColorImageT = RenderColorMesh(faces, vertices, color, cam, ScreenSizeV, ZNearFarV, zoomFactor);
                    
                if type == 2
                    ColorImageT(:,:,1) = flipud(ColorImageT(:,:,1));
                    ColorImageT(:,:,2) = flipud(ColorImageT(:,:,2));
                    ColorImageT(:,:,3) = flipud(ColorImageT(:,:,3));
                end
            end
            
            
            if isDebug
%                 imagesc(DepthImageM); colorbar;
                dImg = mat2gray(DepthImageM);
                if strcmp(colorType, 'color') == 1
                    cImg = double(ColorImageT)./255;
                    imshow([repmat(dImg, 1, 1, 3) cImg]);
                    saveas(h, sprintf('results/syn_test/%0d.png', ii));
                else
                    imagesc(dImg); colorbar;
                end
                pause(0.01);
            else
                
                
                
                
                fwrite(logFile, 1000, 'int64');
                fwrite(logFile, prod(ScreenSizeV)*2, 'int32');
                
                if strcmp(colorType, 'color')==1
                    fwrite(logFile, prod(ScreenSizeV)*3, 'int32');
                else
                    fwrite(logFile, 0, 'int32');
                end
                
                fwrite(logFile, repmat(DepthImageM', 1, 1), 'ushort');
                
                if strcmp(colorType, 'color')==1
                    fwrite(logFile, permute(ColorImageT, [3, 2, 1]), 'int8');
                end
            end
        end
        
        if ~isDebug
            fclose(logFile);
        end

    end
    end

end
