function [ dataSet, poseSet ] = loadEFDatasetBuffer( dataRoot, data_set, bufferDensity)
%LOADEFDATASET Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    bufferDensity = 20;
end

isOverwrite = 0;

poseSet = {};
d_size = length(data_set);

% normVar = [];
% minP = [];

cameraPosMat = {};
cameraPosMat2 = {};

% syn_r = {};
% syn_t = {};

dataSet = cell(d_size, 1);

progressbar2(0);
progressbar2('Loading ply dataset...');
isHaveSubFile = 0;
postFix = sprintf('_d%d', bufferDensity);
for data_idx = 1:d_size
    
%     data_path = sprintf('data/desk/%d', data_idx);
%     data_name = 'MeshedReconstruction';
%     data_ext = 'stl';
%     crop_path = sprintf('%s/crop', data_path);
    
%     fprintf('loading.. %d/%d', data_idx, 1);
    
%     fv{data_idx} = stlread(sprintf('%s/%s_%.2f.%s', crop_path, data_name, dist_set{data_idx}, data_ext));
    
    
    %data3DFilePath = sprintf('%s/%s.klg.ply', dataRoot, data_set{data_idx});
    
    nameParts = strsplit(data_set{data_idx}, '.ply');
    
    subFileName = sprintf('%s%s.mat', nameParts{1}, postFix);
    subFilePath = sprintf('%s/%s', dataRoot, subFileName);
    if exist(subFilePath, 'file') == 2 && isOverwrite == 0
        load(subFilePath);
        isHaveSubFile = 1;
    else
        data3DFilePath = sprintf('%s/%s', dataRoot, data_set{data_idx});
        data3DFilePath
        [tri, pts, data, comments] = ply_read(data3DFilePath, 'tri');    
        isHaveSubFile = 0;
    end
    
    
    
    dataCameraPoseFilePath = sprintf('%s/%s.freiburg', dataRoot, nameParts{1});
    if exist(dataCameraPoseFilePath, 'file') == 2    
        dataCameraPoseFileID = fopen(dataCameraPoseFilePath);
        C = textscan(dataCameraPoseFileID,'%f');
        pData = C{1};
        fclose(dataCameraPoseFileID);
    %     celldisp(C)

        poseSet{data_idx} = [];
        cameraPosMat{data_idx} = [];
        cameraPosMat2{data_idx} = [];
        frameSize = size(pData, 1)/8;
        for fi=0:frameSize-1
            t = pData(fi*8+2:fi*8+4);
    %         pData(fi*8+5:fi*8+8)

            quanternion = pData(fi*8+5:fi*8+8);
            quant_reorder = [quanternion(4) ;quanternion(1:3)];
            r = q_getRotationMatrix(quant_reorder);
            cameraLoc = -r'*t;
            cameraLoc2 = t;

            poseSet{data_idx}(end+1).t = t;
            poseSet{data_idx}(end).r = r;
            poseSet{data_idx}(end).cameraLoc = cameraLoc;
            cameraPosMat{data_idx} = [cameraPosMat{data_idx} cameraLoc];
            cameraPosMat2{data_idx} = [cameraPosMat2{data_idx} cameraLoc2];
        end
    
    end
    
    
    
    
    if isHaveSubFile == 0
        v = [data.vertex.x data.vertex.y data.vertex.z];
        f = tri';
        c = [data.vertex.red data.vertex.green data.vertex.blue];
        n = [data.vertex.nx data.vertex.ny data.vertex.nz];
    
%         fcn_saveUniformSizeModel( v, c, n, dataRoot, nameParts{1}, bufferDensity );
        [ v, c, n, sIndex ] = fcn_uniformSampling( v', c', n', bufferDensity );
        curData = [];
        curData.v = v;
        curData.c = c;
        curData.n = n;
        curData.sIndex = sIndex;
        save(subFilePath, 'curData');
    else
        v = curData.v;
        c = curData.c;
        n = curData.n;
        sIndex = curData.sIndex;
        f = [];
    end
%     r_trans = (rand(3,1)-0.5)/4;
%     r_theta = rand()*pi/16 - pi/32;
%     r_rot = [cos(r_theta) -sin(r_theta) 0 ; sin(r_theta) cos(r_theta) 0 ; 0 0 1];
    
%     syn_t{data_idx} = r_trans;
%     syn_r{data_idx} = r_rot;
    
    %Make intentional error
%     v = bsxfun(@plus, r_rot*v', r_trans)';
    
% 
%     fv{data_idx}.Faces = f;
%     fv{data_idx}.Vertices = v;
%     fv{data_idx}.FaceVertexCData = c./255;
    

    
%     
% %     if data_idx == 1
% %         normVar = sqrt(var(sData')');
% %     end
% %     sData = sData ./ repmat(normVar, 1, size(sData, 2));
%     if data_idx > 1
% %          sData = R*sData;% + repmat(rand(3, 1), 1, size(sData, 2));
% %          sData = sData + repmat(rand(3, 1), 1, size(sData, 2));
%     end
%     
   
    
    
    %Create small version to accerelate occlusion finding
%     fprintf('Generate small version..');
%     fv_small.vertices = fv{data_idx}.Vertices;
%     fv_small.faces = fv{data_idx}.Faces;
%     fv_small = reducepatch(fv_small, 0.01);
%     fprintf('End\n');
%     fv{data_idx}.small = fv_small;
    
    
    %Subsampling
%     [sData, sIndex, sColor] = uniformSubSample(v, 20, c./255);
%     [~, ~, sNormal] = uniformSubSample(v, 20, n);
    
%     V2{end+1} = sData;
%     C2{end+1} = sColor;
%     I2{end+1} = sIndex;
%     Pose_Set{end+1} = cellfun(@(pose) r_rot*pose'+r_trans, poseSet,'uniformoutput',false);
%     PoseMat_Set{end+1} = cellfun(@(pose) [[r_rot*pose(1:3,1:3) ; [0 0 0]] [r_rot*pose(1:end-1,end)+r_trans ;1]], poseMatSet2,'uniformoutput',false);
%     V_all{end+1} = fv{data_idx}.Vertices';
    
    dataSet{data_idx}.v = v';
    dataSet{data_idx}.f = f';
    dataSet{data_idx}.c = (c./255)';
    dataSet{data_idx}.n = (n)';
    dataSet{data_idx}.sIndex = sIndex;
%     dataSet{data_idx}.vs = sData;
%     dataSet{data_idx}.cs = sColor;
%     dataSet{data_idx}.ns = sNormal;
%     dataSet{data_idx}.sIndex = sIndex;
    progressbar2(data_idx / d_size);
end

end

