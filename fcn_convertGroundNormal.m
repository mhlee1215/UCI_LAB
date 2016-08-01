function [] = fcn_convertGroundNormal( dataRoot, inPath, outPath )
%FCN_CONVERTGROUNDNORMAL Summary of this function goes here
%   Detailed explanation goes here

% g = gpuDevice(1);
% reset(g);

disp(sprintf('dataRoot:%s', dataRoot));
disp(sprintf('inPath:%s', inPath));
disp(sprintf('outPath:%s', outPath));

addpath(genpath('/home/mhlee/lab_codes'));
addpath(genpath('/home/mhlee/lab_codes/libs'));
% progressbar2(0);
% maxStep = 5;
% curStep = 1;

data3DFilePath = sprintf('%s/%s', dataRoot, inPath);
data_idx = 1;

disp('Loading ply file...');
% progressbar2(curStep/maxStep); curStep = curStep + 1;

[tri, pts, data0, comments] = ply_read(data3DFilePath, 'tri');
v = [data0.vertex.x data0.vertex.y data0.vertex.z];
f = tri';
c = [data0.vertex.red data0.vertex.green data0.vertex.blue];
n = [data0.vertex.nx data0.vertex.ny data0.vertex.nz];
r = data0.vertex.radius;

dataSet{data_idx}.v = v';
dataSet{data_idx}.f = f';
dataSet{data_idx}.c = (c./255)';
dataSet{data_idx}.n = (n)';

disp('Uniform Sampling...');
% progressbar2(curStep/maxStep); curStep = curStep + 1;

dataSetSampled = loadUniformSampling(dataSet, 10);

for i=1:length(dataSetSampled)
    data.vs{i} = dataSetSampled{i}.v;
    data.ns{i} = dataSetSampled{i}.n;
    data.cs{i} = dataSetSampled{i}.c;
    data.is{i} = dataSetSampled{i}.sIndex;

    data.vb{i} = dataSet{i}.v;
    data.nb{i} = dataSet{i}.n;
    data.cb{i} = dataSet{i}.c;
    data.ib{i} = [];
end

data.vs = data.vs';
data.ns = data.ns';
data.cs = data.cs';
% data.vm = data.vm';
% data.nm = data.nm';
% data.cm = data.cm';
data.vb = data.vb';
data.nb = data.nb';
data.cb = data.cb';



idx = 1;
vForSeg = data.vs{idx}';
cForSeg = data.cs{idx}';
nForSeg = data.ns{idx}';
% iForSeg = data.ns{idx}.sIndex;

% [ segmentedPoints, coloredCloud ] = Segmentation(vForSeg, nForSeg.*0.5, 200, '');
% [segmentIndex, segmentId] = getIndexFromVertices(vForSeg, segmentedPoints);
% 
% Seg.group = segmentIndex;
% Seg.idx = segmentId;
% 
% % Vertex = data.vs{idx};
% % Color = data.cs{idx};
% % SegVertex = coloredCloud(1:3, :);
% % SegColor = coloredCloud(4:6, :);
% 
% segNorm = [];
% for i=1:length(Seg.group)
%     segNorm = [segNorm ; [i sqrt(mean(var(nForSeg(Seg.group{i}, :))))/length(Seg.group{i})]];
% end
% 
% minIdx = 11;
% pclviewer([vForSeg(Seg.group{minIdx},:) cForSeg(Seg.group{minIdx},:)]');
% 
% w = [0.05 1.0 2.0];
% [clustCent,point2cluster,clustMembsCell] = ...
%     MeanShiftCluster([vForSeg.*w(1) cForSeg.*w(2) (nForSeg).*w(3)]', 0.5);
% 
% clustVar = {};
% for i=1:length(clustMembsCell)
%     clustVar{i} = sqrt(mean(var(nForSeg(find(point2cluster == i), :)))) + ...
%         sqrt(mean(var(cForSeg(find(point2cluster == i), :))));
% end
% 
% s=cellfun(@size,clustMembsCell,'uniform',false);
% [trash is]=sortrows(cat(1,s{:}),-[1 2]);
% clusterSorted = clustMembsCell(is);
% % segmentIndex = getIndexFromVertices(V{segId}', clusterSorted);
% % Annotation3DGui(Vertex, Color, SegVertex, SegColor, Seg);
% 
% 
% idxSet = clusterSorted{1};%find(point2cluster == 1);
% pclviewer([vForSeg(idxSet,:) cForSeg(idxSet,:)]');
% 
% 
% pclviewer([vForSeg cForSeg]');

disp('Ransac for fit plane...');
% progressbar2(curStep/maxStep); curStep = curStep + 1;
[B, P, inliers] = ransacfitplane(vForSeg', 0.03);
% pclviewer([vForSeg(inliers,:) cForSeg(inliers,:)]');
meanNormal = -B(1:3);%mean(nForSeg(inliers, :));
R=fcn_RotationFromTwoVectors(meanNormal, [0 0 1]);

Rv = (R*vForSeg')';
zOffset = mean(Rv(inliers, 3));
% pclviewer([(R*vForSeg')' cForSeg]');
% pclviewer([(vForSeg')' cForSeg]');

vForSeg = data.vb{idx}'; 
cForSeg = data.cb{idx}';
nForSeg = data.nb{idx}';

Rv = bsxfun(@minus, (R*vForSeg')' ,[0 0 zOffset]);
Rn = (R*nForSeg')';
d.vertex.x = Rv(:,1);
d.vertex.y = Rv(:,2);
d.vertex.z = Rv(:,3);
d.vertex.red = uint8(cForSeg(:,1).*255);
d.vertex.green = uint8(cForSeg(:,2).*255);
d.vertex.blue = uint8(cForSeg(:,3).*255);
d.vertex.nx = Rn(:,1);
d.vertex.ny = Rn(:,2);
d.vertex.nz = Rn(:,3);
d.vertex.radius = r;

writingPath = sprintf('%s/%s', dataRoot, outPath);
disp(sprintf('Writing Ply...%s', writingPath));
% progressbar2(curStep/maxStep); curStep = curStep + 1;
ply_write(d, writingPath, 'binary_little_endian');


poseFilePath = sprintf('%s/%s.freiburg', dataRoot, inPath);
%     r
     poseMat = fscanf(poseFileID, '%f');
%     poseMat = reshape(poseMat, 17, length(poseMat)/17)';    
%     poseSet = {};
%     poseMatSet = {};
%     for pi = 1:size(poseMat, 1)
%         pM = [poseMat(pi, 1*4+1+1) poseMat(pi, 2*4+1+1) poseMat(pi, 3*4+1+1)];
%         poseMatSet{end+1} = reshape(poseMat(pi, 2:end), 4, 4);
%         poseSet{end+1} = pM;
%     end
%     
%     poseSetMat = cell2mat(poseSet);
%     poseSetMat = reshape(poseSetMat, 3, length(poseSet))';
%     [poseSetMatUnique, ia, ic] = unique(poseSetMat, 'rows');
%     poseSet = mat2cell(poseSetMatUnique, repmat(1, size(poseSetMatUnique,1), 1), 3)';
%     poseMatSet2 = poseMatSet(ia);


disp('Finished.');

end

