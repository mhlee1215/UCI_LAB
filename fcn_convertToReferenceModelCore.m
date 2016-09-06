function [R1, T1] = fcn_convertToReferenceModelCore( points, colors, normals, catName, codeName, fName, isOverwrite )

if nargin < 7
    isOverwrite = 0;
end

R1 = eye(3);
T1 = [0 0 0]';

dataRoot = '/home/mhlee/data_from_odroid';
% catName = 'LAB_1';
refDataRoot = sprintf('%s/match/%s', dataRoot, codeName);
cacheDataRoot = sprintf('%s/cache', dataRoot);
% completeDataRoot = sprintf('%s/complete', dataRoot);
refName = sprintf('%s_ref', catName);
refFeature = sprintf('%s.mat', refName);
refFeaturePath = sprintf('%s/%s', refDataRoot, refFeature);
if exist(refFeaturePath, 'file') ~= 2
    return;
end
load(refFeaturePath, 'feature');
% refFeaturePath = sprintf('%s/%s', refDataRoot, refFeature);


featureNew = {};
params_desc.normalRadius=0.2;
params_desc.searchRadius=1.2;
params_desc.searchK = 0;

v1 = points';
c1 = colors';
n1 = normals';

featureFilePath = sprintf('%s/%s.mat', cacheDataRoot, fName);

if exist(featureFilePath, 'file') == 2 && isOverwrite == 0
     load(featureFilePath);
else
    [k1,~,~] = uniformSubSample(v1, 5, c1);
    descriptors = FPFHSDescriptor(v1, c1, k1, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);

    featureNew.keypoints = k1;
    featureNew.descriptors = descriptors;
    featureNew.v = v1;
    featureNew.c = c1;
    featureNew.n = n1;  
    save(featureFilePath, 'featureNew');
end





% 
% fromIdx = 3;
% toIdx = 1;

v1 = feature.v;
c1 = feature.c;
n1 = feature.n;

v2 = featureNew.v;
n2 = featureNew.n;
c2 = featureNew.c;

% params_desc.normalRadius=0.2;
% params_desc.searchRadius=1.2;
% params_desc.searchK = 0;

% [k1,~,~] = uniformSubSample(v1, 5, c1);
% f1 = FPFHSDescriptor(v1, c1, k1, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);
% mean(sum(descriptors==0,1))
k1 = feature.keypoints;
f1 = feature.descriptors;
 
% [k2,~,~] = uniformSubSample(v2, 5, c2);
% f2 = FPFHSDescriptor(v2, c2, k2, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);
k2 = featureNew.keypoints;
f2 = featureNew.descriptors;

% [ Rf, Tf ] = featureBasedRegistration( v1, c1, n1, v2, c2, n2, k1, f1, k2, f2, 300, true );
% [ Rf, Tf ] = featureBasedRegistration( v1, c1, [], v2, c2, [], k1, f1, k2, f2, 300, true );
[ Rf, Tf ] = featureBasedRegistration( v1, c1, n1, v2, c2, n2, k1, f1, k2, f2, 500, true );
v2t = bsxfun(@plus, Rf*v2', Tf)';
n2t = (Rf*n2')';
% pclviewer([v1 c1 ; v2t c2]');






V = {};
N = {};
C = {};

V{1} = v1';
C{1} = c1';
N{1} = n1';
V{2} = v2t';
C{2} = c2';
N{2} = n2t';


V = V';
C = C';
N = N';

% V = data.vs(idx);
% N = data.ns(idx);
% C = data.cs(idx);

[R,t,X,S,a,pk,T,TAssigned, TXQ] = joint_align(V,N,C, 100);

Rj = cell2mat(T(1,1,end))'*cell2mat(T(2,1,end));
Tj = cell2mat(T(1,1,end))'*cell2mat(T(2,2,end)) - cell2mat(T(1, 2, end));

T1 = Rj * Tf + Tj;
R1 = Rj * Rf;
% 
% v2tt = bsxfun(@plus, R1*v2', T1)';
% n2tt = (R1*n2')';
% pclviewer([v1 c1 ; v2tt c2]');


% params.type = 1;
% params.Assigned = TAssigned;
% params.K = length(TXQ{1, 1, end});
% params.interval = 5;
% % params.marker = marker;
% % params.markerSize = markerSize;
% % params.clrmap = clrmap;
% % params.strIdx = strIdx;
% params.view = [];%[40 54];
% params.pause = 1;
% params.TXQ = TXQ;
% h = figure;
% [TV] = drawTransformation(V, T, params);

% TV = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V,T(:,1,end),T(:,2,end),'uniformoutput',false);

% 
% vaa = [];
% for i = 1:length(TV)
%     vaa = [vaa [ TV{i} ;C{i}]];
% end
% % pclviewer([V{1} V{2};C{1} C{2}]);
% pclviewer(vaa);
% 
% 
% vaa = [];
% 
% vaa = [[feature.v' ; feature.c'] [ bsxfun(@plus, R1*points, T1) ;colors]];
% 
% % pclviewer([V{1} ;C{1}]);
% pclviewer(vaa);
% 




end

