function [] = fcn_modelNormalization( dataRoot, inPath, outPath )


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

disp(sprintf('Loading ply file...%s', data3DFilePath));
tic;
[tri, ~, data0, ~] = ply_read(data3DFilePath, 'tri');
toc;
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
data.vb = data.vb';
data.nb = data.nb';
data.cb = data.cb';



idx = 1;
vForSeg = data.vs{idx}';
cForSeg = data.cs{idx}';
nForSeg = data.ns{idx}';

vForSegBig = data.vb{idx}';
cForSegBig = data.cb{idx}';
nForSegBig = data.nb{idx}';



R = eye(3);
T = [0 0 0]';




disp('Ground Surface normalization..');
[Rground, Tground, inlierIndex] = fcn_convertGroundNormalCore( vForSeg, nForSeg, cForSeg );
% Rv = (Rground*vForSeg')';

Rground
Tground

vGround = bsxfun(@plus, (Rground*vForSeg')' ,Tground');
nGround = (Rground*nForSeg')';

vGroundBig = bsxfun(@plus, (Rground*vForSegBig')' ,Tground');
nGroundBig = (Rground*nForSegBig')';

[ vFloorBig ] = fcn_floorNormalizationCharless( vGroundBig, vGround, inlierIndex );

[vFloor, ~, cFloor] = uniformSubSample(vFloorBig, 10, cForSegBig);
[~, ~, nFloor] = uniformSubSample(vFloorBig, 10, nGroundBig);

% R = R*Rground;
% T = Rground*T + Tground;




% pclviewer([(R*vForSeg')' cForSeg]');
% pclviewer([(vForSeg')' cForSeg]');
% pclviewer([vGround cForSeg]');
%figure; scatter3(vForSeg(:, 1)', vForSeg(:, 2)', vForSeg(:, 3)', 8, cForSeg, 'filled');
% figure; scatter3(vGround(:, 1)', vGround(:, 2)', vGround(:, 3)', 8, cForSeg, 'filled');

fNameParts = strsplit(inPath, '.klg');
fName = fNameParts{1};
nameParts = strsplit(inPath, '-');
catName = nameParts{1};
disp('Convert to fit for reference model');
% refStr = '_norm_ground3';
codeName = 'ref_w_normal';
isOverwrite = 0;
[Rmodel, tmodel] = fcn_convertToReferenceModelCore( vFloor, cFloor, nFloor, catName, codeName, fName, isOverwrite);

Rmodel
tmodel
load(sprintf('/home/mhlee/data_from_odroid/match/%s/%s_ref.mat', codeName, catName));
% vaa = [[feature.v' ; feature.c'] [ bsxfun(@plus, Rmodel*vGround', tmodel) ;cForSeg']];
% pclviewer(vaa);

R = Rmodel*Rground;
T = Rmodel*Tground + tmodel;

% vForSeg = data.vb{idx}'; 
cForSeg = data.cb{idx}';
% nForSeg = data.nb{idx}';

% vForSeg = data.vs{idx}'; 
% cForSeg = data.cs{idx}';
% nForSeg = data.ns{idx}';

vModel = bsxfun(@plus, (Rmodel*vFloorBig')' ,tmodel');
nModel = (Rmodel*nGroundBig')';
Rc = cForSegBig;
% pclviewer([vModel Rc ; feature.v feature.c]')
% pclviewer([vModel Rc]')
% pclviewer([feature.v feature.c]')


% refModel = ply_readWrap(sprintf('/home/mhlee/data_from_odroid/match/%s_ref_t.ply', catName));
% 
% pclviewer([vModel Rc ; refModel.v refModel.c./255]')

d.vertex.x = vModel(:,1);
d.vertex.y = vModel(:,2);
d.vertex.z = vModel(:,3);
d.vertex.red = uint8(cForSeg(:,1).*255);
d.vertex.green = uint8(cForSeg(:,2).*255);
d.vertex.blue = uint8(cForSeg(:,3).*255);
d.vertex.nx = nModel(:,1);
d.vertex.ny = nModel(:,2);
d.vertex.nz = nModel(:,3);
d.vertex.radius = r;

writingPath = sprintf('%s/%s', dataRoot, outPath);
disp(sprintf('Writing Ply...%s', writingPath));
% progressbar2(curStep/maxStep); curStep = curStep + 1;
ply_write(d, writingPath, 'binary_little_endian');


normInfoFilePath = sprintf('%s/%s.info', dataRoot, outPath);
disp(sprintf('Writing Info...%s', normInfoFilePath));
normInfoFile = fopen(normInfoFilePath, 'w');
q = q_getFromRotationMatrix(R);
fprintf(normInfoFile, '%f %f %f %f %f', q(1), q(2), q(3), q(4), T);
fclose(normInfoFile);
    
disp('Finished.');

end

