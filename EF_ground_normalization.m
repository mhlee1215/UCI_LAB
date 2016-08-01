
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
iForSeg = data.ns{idx}.sIndex;

[ segmentedPoints, coloredCloud ] = Segmentation(vForSeg, nForSeg.*0.5, 200, '');
[segmentIndex, segmentId] = getIndexFromVertices(vForSeg, segmentedPoints);

Seg.group = segmentIndex;
Seg.idx = segmentId;

% Vertex = data.vs{idx};
% Color = data.cs{idx};
% SegVertex = coloredCloud(1:3, :);
% SegColor = coloredCloud(4:6, :);

segNorm = [];
for i=1:length(Seg.group)
    segNorm = [segNorm ; [i sqrt(mean(var(nForSeg(Seg.group{i}, :))))/length(Seg.group{i})]];
end

minIdx = 11;
pclviewer([vForSeg(Seg.group{minIdx},:) cForSeg(Seg.group{minIdx},:)]');

w = [0.05 1.0 2.0];
[clustCent,point2cluster,clustMembsCell] = ...
        MeanShiftCluster([vForSeg.*w(1) cForSeg.*w(2) (nForSeg).*w(3)]', 0.5);

clustVar = {};
for i=1:length(clustMembsCell)
    clustVar{i} = sqrt(mean(var(nForSeg(find(point2cluster == i), :)))) + ...
        sqrt(mean(var(cForSeg(find(point2cluster == i), :))));
end
    
s=cellfun(@size,clustMembsCell,'uniform',false);
[trash is]=sortrows(cat(1,s{:}),-[1 2]);
clusterSorted = clustMembsCell(is);
% segmentIndex = getIndexFromVertices(V{segId}', clusterSorted);    
% Annotation3DGui(Vertex, Color, SegVertex, SegColor, Seg);


idxSet = clusterSorted{1};%find(point2cluster == 1);
pclviewer([vForSeg(idxSet,:) cForSeg(idxSet,:)]');


pclviewer([vForSeg cForSeg]');
[B, P, inliers] = ransacfitplane(vForSeg', 0.03);
pclviewer([vForSeg(inliers,:) cForSeg(inliers,:)]');
meanNormal = mean(nForSeg(inliers, :));
R=fcn_RotationFromTwoVectors(meanNormal, [0 0 1]);

pclviewer([(R*vForSeg')' cForSeg]');
pclviewer([(vForSeg')' cForSeg]');

Rv = (R*vForSeg')';
Rn = (R*nForSeg')';
d.vertex.x = Rv(:,1);
d.vertex.y = Rv(:,2);
d.vertex.z = Rv(:,3);
d.vertex.red = uint8(cForSeg(:,1).*255);
d.vertex.green = uint8(cForSeg(:,2).*255);
d.vertex.blue = uint8(cForSeg(:,3).*255);
d.vertex.nx = Rn(:,1);
d.vertex.ny = Rn(:,1);
d.vertex.nz = Rn(:,1);
ply_write(d, 'testSavePly.ply');




