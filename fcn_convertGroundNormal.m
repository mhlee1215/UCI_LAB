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

[Rground] = fcn_convertGroundNormalCore( vForSeg, nForSeg );

Rv = (Rground*vForSeg')';
zOffset = mean(Rv(inliers, 3));
% pclviewer([(R*vForSeg')' cForSeg]');
% pclviewer([(vForSeg')' cForSeg]');
%figure; scatter3(vForSeg(:, 1)', vForSeg(:, 2)', vForSeg(:, 3)', 8, cForSeg, 'filled');
% figure; scatter3(Rv(:, 1)', Rv(:, 2)', Rv(:, 3)', 8, cForSeg, 'filled');



vForSeg = data.vb{idx}'; 
cForSeg = data.cb{idx}';
nForSeg = data.nb{idx}';

Rv = bsxfun(@minus, (Rground*vForSeg')' ,[0 0 zOffset]);
Rn = (Rground*nForSeg')';
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


normInfoFilePath = sprintf('%s/%s.info', dataRoot, outPath);
disp(sprintf('Writing Info...%s', normInfoFilePath));
normInfoFile = fopen(normInfoFilePath, 'w');
q = q_getFromRotationMatrix(Rground);
fprintf(normInfoFile, '%f %f %f %f %f', q(1), q(2), q(3), q(4), zOffset);
fclose(normInfoFile);
    
disp('Finished.');

end

