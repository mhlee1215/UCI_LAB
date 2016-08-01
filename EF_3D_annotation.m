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
dataRoot = '/home/mhlee/data_from_odroid/complete';
% data_set = {'2016-02-09.01','2016-02-09.02','2016-02-09.03'};
% data_set = {'2016-04-29.03', '2016-04-29.05', '2016-04-29.06', '2016-04-29.00','2016-04-29.08'};
%Wrong init pos : '2016-04-29.06'
% data_set = {'2016-04-29.00', '2016-04-29.01', '2016-04-29.02', '2016-04-29.03', '2016-04-29.05', '2016-04-29.08'};%, '2016-04-29.05'};%, '2016-04-29.00','2016-04-29.08'};
data_set = {'MAC_TEST-2016-07-17_23_11'};
d_size = length(data_set);
% 
% 
[dataSet, poseSet] = loadEFDataset( dataRoot, data_set);
dataSetSampled = loadUniformSampling(dataSet, 10);
dataSetSampled2 = loadUniformSampling(dataSet, 20);
% 


% V = {};
% N = {};
% C = {};
% Vall = {};
% Nall = {};
% Call = {};
for i=1:length(dataSetSampled)
    data.vs{i} = dataSetSampled{i}.v;
    data.ns{i} = dataSetSampled{i}.n;
    data.cs{i} = dataSetSampled{i}.c;
    data.is{i} = dataSetSampled{i}.sIndex;
    data.vm{i} = dataSetSampled2{i}.v;
    data.nm{i} = dataSetSampled2{i}.n;
    data.cm{i} = dataSetSampled2{i}.c;
    data.im{i} = dataSetSampled2{i}.sIndex;
    data.vb{i} = dataSet{i}.v;
    data.nb{i} = dataSet{i}.n;
    data.cb{i} = dataSet{i}.c;
    data.ib{i} = [];
end

data.vs = data.vs';
data.ns = data.ns';
data.cs = data.cs';
data.vm = data.vm';
data.nm = data.nm';
data.cm = data.cm';
data.vb = data.vb';
data.nb = data.nb';
data.cb = data.cb';



idx = 1;
vForSeg = data.vs{idx}';
cForSeg = data.cs{idx}';
nForSeg = data.ns{idx}';
iForSeg = data.is{idx};


[ segmentedPoints, coloredCloud ] = Segmentation(vForSeg, cForSeg, 200, '');
[segmentIndex, segmentId] = getIndexFromVertices(vForSeg, segmentedPoints');

Seg.group = segmentIndex;
Seg.idx = segmentId;

AnnotParams.Vertex = vForSeg';
AnnotParams.Color = cForSeg';
AnnotParams.SegVertex = coloredCloud(1:3, :);
AnnotParams.SegColor = coloredCloud(4:6, :);
AnnotParams.Seg = Seg;
% Vertex = vForSeg';
% Color = cForSeg';
% SegVertex = coloredCloud(1:3, :);
% SegColor = coloredCloud(4:6, :);


Annotation3DGui(AnnotParams);




