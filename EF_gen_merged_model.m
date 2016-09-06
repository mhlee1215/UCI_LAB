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



for mi = 1:length(category_set)

modelName = category_set{mi};%'LAB_1';
fileList = dir(sprintf('%s%s-*cvt.ply', dataRoot, modelName));

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
% dataSetSampled2 = loadUniformSampling(dataSet, 20);
% 


% V = {};
% N = {};
% C = {};
% Vall = {};
% Nall = {};
% Call = {};
data = {};
% MarkColor = {};
for i=1:length(dataSet)

    data.vb{i} = dataSet{i}.v;
    data.nb{i} = dataSet{i}.n;
    data.cb{i} = dataSet{i}.c;
%     data.ib{i} = dataSet{i}.sIndex;
    
    data.vs{i} = dataSetSmall{i}.v;
    data.ns{i} = dataSetSmall{i}.n;
    data.cs{i} = dataSetSmall{i}.c;
    data.is{i} = dataSetSmall{i}.sIndex;
end

% data.vs = data.vs';
% data.ns = data.ns';
% data.cs = data.cs';
% data.vm = data.vm';
% data.nm = data.nm';
% data.cm = data.cm';
data.vb = data.vb';
data.nb = data.nb';
data.cb = data.cb';
data.vs = data.vs';
data.ns = data.ns';
data.cs = data.cs';

V = data.vs;
C = data.cs;
N = data.ns;

[R,t,X,S,a,pk,T,TAssigned, TXQ, vis] = joint_align(V, N, C, 50, 2);

TV = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V,T(:,1,end),T(:,2,end),'uniformoutput',false);
TVb = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),data.vb,T(:,1,end),T(:,2,end),'uniformoutput',false);

visMat = cell2mat(vis');
visMat = visMat > 0.06;

visMatCol = ones(1, size(visMat, 2));
for i = 1:size(visMat, 1)
    visMatCol = visMatCol .* visMat(i, :);
end

visMatCol = visMatCol + 1;


SEG_color = [[0.3 0.3 0.3];[0.9 0.1 0.1]; [0.1 0.9 0.1] ; [0.1 0.1 0.9]];

pcl_model = [];
Rv = [];
Rc = [];
for i=1:length(V)
    clustAssin = cell2mat(TAssigned(i,end));
    dynamicObjVec = visMatCol(clustAssin);    
    pcl_model = [pcl_model [TV{i} ;SEG_color(dynamicObjVec, :)']];
    
    Rv = [Rv TV{i}];
    Rc = [Rc SEG_color(dynamicObjVec, :)'];
end
% pclviewer(pcl_model);
fileName = sprintf('%s_merged_static2', modelName);
dataRoot2 = '/home/mhlee/data_from_odroid/merged/';
fcn_saveUniformSizeModel( Rv', Rc', [], dataRoot2, fileName, 0 );

pcl_model = [];
Rv = [];
Rc = [];
for i=1:length(V)
    clustAssin = cell2mat(TAssigned(i,end));
    dynamicObjVec = visMatCol(clustAssin);    
    pcl_model = [pcl_model [TVb{i} ;data.cb{i}]];
    Rv = [Rv TVb{i}];
    Rc = [Rc data.cb{i}];
end
% pclviewer(pcl_model);
fileName = sprintf('%s_merged_rgb2', modelName);
dataRoot2 = '/home/mhlee/data_from_odroid/merged/';
fcn_saveUniformSizeModel( Rv', Rc', [], dataRoot2, fileName, 0 );

pcl_model = [];
Rv = [];
Rc = [];
for i=1:length(V)
    pcl_model = [pcl_model [TVb{i} ; repmat(MarkColor{i}', 1, size(TVb{i}, 2))]];
    Rv = [Rv TVb{i}];
    Rc = [Rc repmat(MarkColor{i}', 1, size(TVb{i}, 2))];
end
% pclviewer(pcl_model);
fileName = sprintf('%s_merged_color2', modelName);
dataRoot2 = '/home/mhlee/data_from_odroid/merged/';
fcn_saveUniformSizeModel( Rv', Rc', [], dataRoot2, fileName, 0 );
end

pcl_model = [];
for i=1:length(V)
    pcl_model = [pcl_model [TVb{i} ; repmat(MarkColor{i}', 1, size(data.vb{i}, 2))]];
end
pclviewer(pcl_model);
% 
% 
% pclviewer([TVb{5} ; data.cb{5}]);