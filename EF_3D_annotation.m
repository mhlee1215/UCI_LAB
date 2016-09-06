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
% run('libs2/vlfeat-0.9.20-bin/vlfeat-0.9.20/toolbox/vl_setup');


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
dataRoot = '/home/mhlee/data_from_odroid/complete/';




%Read All category set
fileListAll = dir(sprintf('%sLAB*cvt.ply', dataRoot));
category_set = {};
strCategory = '';
for i=1:length(fileListAll)
    nameParts = strsplit(fileListAll(i).name, '-');
    %     nameParts2 = strsplit(nameParts{4}, '_');
    curStrCategory = nameParts{1};%sprintf('%s-%s-%s', , nameParts{3}, nameParts2{1});
    if strcmp(strCategory, curStrCategory) == 1
        continue;
    end
    %     data_set{end+1} = fileList(i).name;
    category_set{end+1} = curStrCategory;
    strCategory = curStrCategory;
end

colors = distinguishable_colors(30);
MarkColor = {};
for i=1:30
    MarkColor{i} = colors(i,:);
end

% for mi = 1:length(category_set)
mi=1;
modelName = category_set{mi};%'LAB_1';
fileList = dir(sprintf('%s%s*cvt.ply', dataRoot, modelName));

date_set = {};
data_set = {};
strDate = '';
for i=1:length(fileList)
    nameParts = strsplit(fileList(i).name, '-');
    nameParts2 = strsplit(nameParts{4}, '_');
    curStrDate = sprintf('%s-%s-%s', nameParts{2}, nameParts{3}, nameParts2{1});
    if strcmp(strDate, curStrDate) == 1
        continue;
    end
    data_set{end+1} = fileList(i).name;
    date_set{end+1} = curStrDate;
    strDate = curStrDate;
end
[dataSet, poseSet] = loadEFDataset( dataRoot, data_set);
% [dataSetSmall, ~] = loadEFDatasetBuffer( dataRoot, data_set, 10);
dataSetSmall = loadUniformSampling(dataSet, 10);
dataSetSampled2 = loadUniformSampling(dataSet, 20);
%


plyData = containers.Map('KeyType', 'int32', 'ValueType', 'any');
for density = [0, 10, 20]
   
    
    if density == 0
        d = dataSet;
    else
        d = loadUniformSampling(dataSet, density);
    end
    
    data = {};
    for i=1:length(d)
        data.v{i} = d{i}.v;
        data.n{i} = d{i}.n;
        data.c{i} = d{i}.c;
        
        if isfield(d{i}, 'sIndex')
            data.si{i} = d{i}.sIndex;
        end
    end
    
    data.v = data.v';
    data.n = data.n';
    data.c = data.c';
    
    if isfield(data, 'si')
        data.si = data.si';
    end
    
    plyData(density) = data;
end


curData = plyData(10);
V = curData.v;
C = curData.c;
N = curData.n;

[R,t,X,S,a,pk,T,TAssigned, TXQ, vis] = joint_align(V, N, C, 10, 5);

% TV = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V,T(:,1,end),T(:,2,end),'uniformoutput',false);
% TVb = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),data.vb,T(:,1,end),T(:,2,end),'uniformoutput',false);

visMat = cell2mat(vis');
visMat = visMat > 0.06;

visMatCol = ones(1, size(visMat, 2));
for i = 1:size(visMat, 1)
    visMatCol = visMatCol .* visMat(i, :);
end

visMatCol = visMatCol + 1;



segData = containers.Map('KeyType', 'int32', 'ValueType', 'any');
segDensity = 20;
curData = plyData(segDensity);
for idx = 1:length(curData.v)
    V = curData.v;
    TV = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V,T(:,1,end),T(:,2,end),'uniformoutput',false);
    vForSeg = TV{idx}';%data.vs{idx}';
    cForSeg = curData.c{idx}';
    nForSeg = curData.n{idx}';
%     iForSeg = data.is{idx}';


    [ segmentedPoints, coloredCloud ] = Segmentation(vForSeg, cForSeg, 200, '');
    [segmentIndex, segmentId] = getIndexFromVertices(vForSeg, segmentedPoints');

    Seg{idx}.group = segmentIndex;
    Seg{idx}.idx = segmentId;
    Seg{idx}.SegVertex = coloredCloud(1:3, :);
    Seg{idx}.SegColor = coloredCloud(4:6, :);
    
end
segData(segDensity) = Seg;

AnnotParams.T = T;
AnnotParams.plyData = plyData;
AnnotParams.Seg = segData;
AnnotParams.TAssigned = TAssigned;
AnnotParams.date_set = date_set;
AnnotParams.name_set = data_set;

Annotation3DGui(AnnotParams);




