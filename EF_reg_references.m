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
colors = distinguishable_colors(30, [0.3 0.3 0.3]);
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
originName = 'LAB_1_ref.ply'; %Which is never move
netInfo = containers.Map();
netInfo('LAB_3_ref.ply') = 'LAB_2_ref.ply';
netInfo('LAB_2_ref.ply') = 'LAB_1_ref.ply';
% netInfo('LAB_3_ref.ply') = 'LAB_3a_ref.ply';
netInfo('LAB_3a_ref.ply') = 'LAB_3_ref.ply';
netInfo('LAB_3b_ref.ply') = 'LAB_1_ref.ply';
netInfo('LAB_4a_ref.ply') = 'LAB_4_ref.ply';
netInfo('LAB_4_ref.ply') = 'LAB_2_ref.ply';

nameToData = containers.Map();
nameToId = containers.Map();
idToData = containers.Map('KeyType', 'int32', 'ValueType', 'any');

codeName = 'ref_w_normal';
dataRoot = '/home/mhlee/data_from_odroid/match';
dataOutRoot = sprintf('%s/%s', dataRoot, codeName);
mkdir(dataOutRoot);
matchInfo = customXml2struct('/home/mhlee/data_from_odroid/match/match_22_ref.mlp');
data_set = {};
data_set.name = {};
data_set.R = {};
data_set.t = {};
for i=1:length(matchInfo.MeshLabProject.MeshGroup.MLMesh)
    data_set.name{i} = matchInfo.MeshLabProject.MeshGroup.MLMesh{i}.Attributes.filename;
    MatText = matchInfo.MeshLabProject.MeshGroup.MLMesh{i}.MLMatrix44.Text;
    
    curR = [];
    MatRows = strsplit(strtrim(MatText), '\n');
    for j = 1:length(MatRows)
        elements = strsplit(strtrim(MatRows{j}), ' ');
        for k = 1:length(elements)
            curR(end+1) = str2num(elements{k});
        end
    end
    curR = reshape(curR, 4, 4)';
    data_set.R{i} = curR(1:3, 1:3);
    data_set.t{i} = curR(1:3, 4);
end

% data_set = {'LAB_1-2016-07-22_12_32.klg_cvt', 'LAB_2-2016-07-22_12_41.klg_cvt'};

% d_size = length(data_set);
% 
% 
[dataSet, poseSet] = loadEFDataset( dataRoot, data_set.name);


for i=1:length(dataSet)
    smallDotIdx = 1:size(dataSet{i}.v, 2);%find(dataSet{i}.radius<0.008);
    dataSet{i}.v = dataSet{i}.v(:, smallDotIdx);
    dataSet{i}.c = dataSet{i}.c(:, smallDotIdx);
    dataSet{i}.n = dataSet{i}.n(:, smallDotIdx);
    dataSet{i}.radius = dataSet{i}.radius(smallDotIdx);
end

dataSetSampled = loadUniformSampling(dataSet, 10);

% V = {};
% N = {};
% C = {};
% Vall = {};
% Nall = {};
% Call = {};
for i=1:length(dataSetSampled)
    
    R0 = data_set.R{i};
    T0 = data_set.t{i};
    R1 = eye(3);
    T1 = [0 0 0]';
    
    
    %Apply manually generated R and T
    curData.vs = bsxfun(@plus, R0*dataSetSampled{i}.v, T0);
    curData.ns = R0*dataSetSampled{i}.n;
    curData.cs = dataSetSampled{i}.c;
    curData.is = dataSetSampled{i}.sIndex;
    
    %Apply ground surface normalization
    [Rground, Tground, inlierIndex] = fcn_convertGroundNormalCore( curData.vs', curData.ns', curData.cs' );
    
    R1 = Rground * R0;
    T1 = Rground * T0 + Tground;
    
    vGround = bsxfun(@plus, (Rground*curData.vs)' ,Tground');
    nGround = (Rground*curData.ns)';
    
    %Apply surface rectification (= all zero z coordination)
    [ vFloor ] = fcn_floorNormalization( vGround, vGround, inlierIndex );
    [ vFloor2 ] = fcn_floorNormalizationCharless( vGround, vGround, inlierIndex );
    
    curData.vsf = vFloor';
    curData.vsf2 = vFloor2';
    curData.nsf = nGround';
    curData.vsg = vGround';
    curData.nsg = nGround';
    
%     pclviewer([curData.vs(:, :) ;curData.cs(:, :)]);
%     pclviewer([curData.vs(:, inlierIndex) ;curData.cs(:, inlierIndex)]);
    
%     vv = [];
%     for i=1:10000
%         r = randperm(length(inlierIndex));
%         v = mean(std(curData.cs(:, r(1:3))'));
%         vv = [vv v];
%     end
%     figure; histogram(vv);
        %     pclviewer([curData.vs ;curData.cs]);
%     figure; scatter3(vGround(:, 1)', vGround(:, 2)', vGround(:, 3)', 8, curData.cs', 'filled');
%     figure; scatter3(vFloor(:, 1)', vFloor(:, 2)', vFloor(:, 3)', 8, curData.cs', 'filled');    
    
    
    curData.vb = bsxfun(@plus, R0*dataSet{i}.v, T0);
    curData.nb = R0*dataSet{i}.n;
    vGroundBig = bsxfun(@plus, (Rground*curData.vb)' ,Tground');
    nGroundBig = (Rground*curData.nb)';
    [ vFloorBig ] = fcn_floorNormalization( vGroundBig, vGround, inlierIndex );
    [ vFloorBig2 ] = fcn_floorNormalizationCharless( vGroundBig, vGround, inlierIndex );
    
    curData.vbf = vFloorBig';
    curData.vbf2 = vFloorBig2';
    curData.nbf = nGroundBig';   
    
    curData.vbg = vGroundBig';   
    curData.nbg = nGroundBig';   
    
    curData.cb = dataSet{i}.c;
    curData.ib = [];
    curData.radius = dataSet{i}.radius;
    
%     pclviewer([curData.vb ;curData.cb]);
    
    nameToData(data_set.name{i}) = curData;
    idToData(i) = curData;
    nameToId(data_set.name{i}) = i;
end

% pclviewer([idToData(7).vbf ;idToData(7).cb]);
% pclviewer([idToData(7).vb ;idToData(7).cb]);
% pclviewer([idToData(5).vb ;idToData(5).cb]);


% idx = 1;
% pclviewer([idToData(idx).vbf' repmat([1 0 0], size(idToData(idx).vbf', 1), 1) ...
%     ; idToData(idx).vbg' repmat([0 1 0], size(idToData(idx).vbg', 1), 1) ...
%     ; idToData(idx).vbf2' repmat([0 1 1], size(idToData(idx).vbg', 1), 1)]')
% 
% for idx = 1:7
%     fileName = sprintf('warp_diff_%d', idx);
%     CC = [repmat([1 0 0], size(idToData(idx).vbf', 1), 1) ; repmat([0 1 0], size(idToData(idx).vbg', 1), 1) ; repmat([0 1 1], size(idToData(idx).vbf2', 1), 1)];
%     fcn_saveUniformSizeModel( [idToData(idx).vbf'; idToData(idx).vbg' ; idToData(idx).vbf2'], CC, [idToData(idx).nbf'; idToData(idx).nbg' ; idToData(idx).nbf'], dataRoot, fileName, 50 );
% end

% vaa = [];
% for i = 1:3
%     vaa = [vaa [ idToData(i).vs ;idToData(i).cs]];
% end
% % pclviewer([V{1} ;C{1}]);
% pclviewer(vaa);
% 
% 
% return;
% save('./cache/EF_reg_references_pre.mat', '*', '-v7.3');

transMatToOrigin = containers.Map();
for i=1:length(dataSetSampled)
    if strcmp(data_set.name{i}, originName) == 1
        trans.R = eye(3);
        trans.T = [0 0 0]';
        transMatToOrigin(data_set.name{i}) = trans;
        continue;
    end
    
%     if  strcmp(data_set.name{i}, 'LAB_4a_ref.ply') == 1 &&
%         strcmp(data_set.name{i}, data_set.name{i}) == 1
%         
%     end
    
    curModel = data_set.name{i};
    if transMatToOrigin.isKey(curModel)
        R1 = transMatToOrigin(curModel).R;
        T1 = transMatToOrigin(curModel).T;
    else
        R1 = eye(3);
        T1 = [0 0 0]';
        while strcmp(curModel, originName) ~= 1
            
            if transMatToOrigin.isKey(curModel)
                R0 = transMatToOrigin(curModel).R;
                T0 = transMatToOrigin(curModel).T;
                R1 = R0 * R1;
                T1 = R0 * T1 + T0;
                break;
            else
            
                fprintf(sprintf('%s -> %s, ', curModel, netInfo(curModel)));
                id2 = nameToId(curModel);
                id1 = nameToId(netInfo(curModel));

                idx = [id1, id2];

                V = {};
                N = {};
                C = {};
                Vb = {};
                Nb = {};
                Cb = {};
                CbSeg = {};
                Rb = {};    
                IS = {};

                for ii=idx
                    V{end+1} = idToData(ii).vsf2;
                    C{end+1} = idToData(ii).cs;
                    N{end+1} = idToData(ii).nsf;
                    IS{end+1} = idToData(ii).is;
                    Vb{end+1} = idToData(ii).vbf2;
                    Cb{end+1} = idToData(ii).cb;
                    CbSeg{end+1} = repmat(colors(ii, :)', 1, size(idToData(ii).cb, 2));
                    Nb{end+1} = idToData(ii).nbf;
                    Rb{end+1} = idToData(ii).radius;
                end

                V = V';
                C = C';
                N = N';
                Vb = Vb';
                Cb = Cb';
                Nb = Nb';
                Rb = Rb';
                IS = IS';

                [R,t,X,S,a,pk,T,TAssigned, TXQ, vis] = joint_align(V,N,C, 50, 5);

                
                %TTT = cell2mat(T);
                %figure; plot(squeeze(TTT(1, :, :))')
                
%                 params.type = 1;
%                 params.Assigned = TAssigned;
%                 params.K = length(TXQ{1, 1, end});
%                 params.interval = 5;
%                 params.view = [];%[40 54];
%                 params.pause = 1;
%                 params.TXQ = TXQ;
%                 h = figure;
%                 drawTransformation(V, T, params);
                
%                 TV = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V,T(:,1,end),T(:,2,end),'uniformoutput',false);
%                 TN = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),N,T(:,1,end),T(:,2,end),'uniformoutput',false);
%                 TVb = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),Vb,T(:,1,end),T(:,2,end),'uniformoutput',false);
%                 TNb = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),Nb,T(:,1,end),T(:,2,end),'uniformoutput',false);
%                 
%                 
%                 pcl_model = [];
%                 for i=1:length(V)
%                     pcl_model = [pcl_model [TVb{i}(:, :) ; Cb{i}]];
% %                     pcl_model = [pcl_model [TVb{i}(:, :) ; repmat(colors(i,:)', 1, size(TVb{i}, 2))]];
% %                     RcSeg = [RcSeg ; repmat(colors(i,:), size(Vt, 2), 1)];
%                 end
%                 pclviewer(pcl_model);
%                 fileName = 'flat_with_normal_RGB';
%                 fcn_saveUniformSizeModel( TVb, Cb, Nb, dataRoot, fileName, 0);
%                 fileName = 'flat_with_normal_Segment';
%                 fcn_saveUniformSizeModel( TVb, CbSeg, Nb, dataRoot, fileName, 0);
                
%                 pcl_model2 = pcl_model;
%                 pclviewer(pcl_model2);
                
%                 visMat = cell2mat(vis');
%                 visMat = visMat > 0.06;
% 
%                 visMatCol = ones(1, size(visMat, 2));
%                 for i = 1:size(visMat, 1)
%                     visMatCol = visMatCol .* visMat(i, :);
%                 end
%                 visMatCol = visMatCol + 1;
                
                
%                 
%                 
%                 
%                 
%                 
%                 V2 = {};
%                 C2 = {};
%                 N2 = {};
% %               
%                 rMax = -[999999 999999 999999];
%                 rMin = -rMax;
%                 for i=1:length(V)
%                     clustAssin = cell2mat(TAssigned(i,end));
%                     dynamicObjVec = visMatCol(clustAssin);    
%                     
%                     overlapIdx = find(dynamicObjVec == 2);
%                     if i==1
%                         rMax = max(TV{i}(:, overlapIdx)');
%                         rMin = min(TV{i}(:, overlapIdx)');
%                     else
%                         rMax = min(rMax, max(TV{i}(:, overlapIdx)'));
%                         rMin = max(rMin, min(TV{i}(:, overlapIdx)'));    
%                     end
%                     
%                 end
%                 
%                 V2 = {};
%                 C2 = {};
%                 N2 = {};
%                 pcl_model = [];
%                 for i=1:length(TV)
%                     VTarget = TVb;
%                     insideIdx = (VTarget{i}(1,:) > rMin(1)) .* (VTarget{i}(2,:) > rMin(2)) .* (VTarget{i}(3,:) > rMin(3));
%                     insideIdx = insideIdx.* (VTarget{i}(1,:) < rMax(1)) .* (VTarget{i}(2,:) < rMax(2)) .* (VTarget{i}(3,:) < rMax(3));
%                     insideIdx = find(insideIdx);
%                     insideIdx = insideIdx(1:5:end);
%                     V2{i} = VTarget{i}(:, insideIdx);
%                     C2{i} = Cb{i}(:, insideIdx);
%                     N2{i} = TNb{i}(:, insideIdx);
%                     
%                     pcl_model = [pcl_model [V2{i}(:, :) ; C2{i}(:, :)]];
%                     pcl_model = [pcl_model [V2{i}(:, :) ; repmat(colors(i, :)', 1, size(V2{i}, 2))]];
%                      
%                 end
%                 pclviewer(pcl_model);
%                 V2 = V2';
%                 C2 = C2';
%                 N2 = N2';
%                 
%                 
%                 
%                 
%                 
%                 
%                 res = 0.3;
%                 [xx,yy,zz] = meshgrid(rMin(1):res:rMax(1),rMin(2):res:rMax(2), rMin(3):res:rMax(3));
%                 v_seed = [xx(:), yy(:), zz(:)]';
%                 v_index = zeros(length(xx(:)), 1);
%                 
%                 kdtree_single = vl_kdtreebuild(v_seed);
%                 [index_single, distance_single] = vl_kdtreequery(kdtree_single, v_seed, V2{1}, 'NumNeighbors', 1, 'MaxComparisons', 10);
%                 
%                 validIdx = find(distance_single < res);
%                 v_index(index_single(validIdx)) = v_index(index_single(validIdx)) + 1;
%                 
% %                 kdtree_single = vl_kdtreebuild(v_seed);
%                 [index_single2, distance_single2] = vl_kdtreequery(kdtree_single, v_seed, V2{2}, 'NumNeighbors', 1, 'MaxComparisons', 10);
%                 
%                 validIdx2 = find(distance_single2 < res);
%                 v_index(index_single2(validIdx2)) = v_index(index_single2(validIdx2)) + 1;
%                 
%                 V3 = {};
%                 C3 = {};
%                 N3 = {};
%                 pcl_model2 = [];
% %                 for i=1:length(TV)
%                     Lia1 = ismember(index_single, find(v_index == 2));
%                     V3{1} = V2{1}(:, find(Lia1));
%                     C3{1} = C2{1}(:, find(Lia1));
%                     N3{1} = N2{1}(:, find(Lia1));
%                     pcl_model2 = [pcl_model2 [V3{1}(:, :) ; repmat(colors(1, :)', 1, length(find(Lia1)))]];
%                     Lia2 = ismember(index_single2, find(v_index == 2));
%                     V3{2} = bsxfun(@plus, V2{2}(:, find(Lia2)), [0.5 0.5 0]');
%                     C3{2} = C2{2}(:, find(Lia2));
%                     N3{2} = N2{2}(:, find(Lia2));
%                     pcl_model2 = [pcl_model2 [V3{2}(:, :) ; repmat(colors(2, :)', 1, length(find(Lia2)))]];
% %                 end
% %                 pclviewer(pcl_model2);
% %                 
%                 V3 = V3';
%                 C3 = C3';
%                 N3 = N3';
%                 
%                 
%                 
%                 
%                 [R2,t2,X2,S2,a2,pk2,T2,TAssigned2, TXQ2, vis2] = joint_align(V3,N3,C3, 300, 8);
% %                 [R2,t2,X2,S2,a2,pk2,T2,TAssigned2, TXQ2, vis2] = joint_align(V2,N2,C2, 150, 10);
% 
%                 
%                 TV2 = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V2,T2(:,1,end),T2(:,2,end),'uniformoutput',false);
% %                 TVb2 = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),Vb,T2(:,1,end),T2(:,2,end),'uniformoutput',false);
%                 
%                 params.type = 1;
%                 params.Assigned = TAssigned2;
%                 params.K = length(TXQ2{1, 1, end});
%                 params.interval = 3;
%                 params.view = [];%[40 54];
%                 params.pause = 0.1;
%                 params.TXQ = TXQ2;
%                 h = figure;
%                 drawTransformation(V2, T2, params);
%                 
%                 
%                 pcl_model = [];
%                 for i=1:length(V)
%                     pcl_model = [pcl_model [TVb{i}(:, :) ; Cb{i}]];
% 
%                 end
%                 pclviewer(pcl_model);
%                 SEG_color = [[0.3 0.3 0.3];[0.9 0.1 0.1]; [0.1 0.9 0.1] ; [0.1 0.1 0.9]];
% 
%                 pcl_model = [];
%                 Rv = [];
%                 Rc = [];
%                 
%                 for i=1:length(V)
%                     clustAssin = cell2mat(TAssigned(i,end));
%                     dynamicObjVec = visMatCol(clustAssin);    
% %                     pcl_model = [pcl_model [TV{i} ;SEG_color(dynamicObjVec, :)']];
%                     smallDotIdx = find(Rb{i}<0.008);
%                     radius = Rb{i};
%                     radius = radius ./ max(radius);
% %                     pcl_model = [pcl_model [TVb{i}(:, smallDotIdx) ;SEG_color(dynamicObjVec(IS{i}(smallDotIdx)), :)']];
%                     rd = radius(smallDotIdx);
%                     pcl_model = [pcl_model [TVb{i}(:, smallDotIdx) ; [rd' ; rd' ;rd']]];
% 
%                     Rv = [Rv TV{i}];
%                     Rc = [Rc SEG_color(dynamicObjVec, :)'];
%                 end
%                 pclviewer(pcl_model);
%                 
%                 pcl_model = [];
%                 Rv = [];
%                 Rc = [];
%                 for i=1:length(V)
%                     pcl_model = [pcl_model [TVb{i} ; repmat(colors(i, :)', 1, size(TVb{i}, 2))]];
% %                     Rv = [Rv TVb{i}];
% %                     Rc = [Rc repmat(colors(:, i), 1, size(TVb{i}, 2))];
%                 end
%                 pclviewer(pcl_model);
%                 
%                 pcl_model = [];
%                 Rv = [];
%                 Rc = [];
%                 for i=1:length(V)
%                     pcl_model = [pcl_model [TVb2{i} ; repmat(colors(i, :)', 1, size(TVb2{i}, 2))]];
% %                     Rv = [Rv TVb{i}];
% %                     Rc = [Rc repmat(colors(:, i), 1, size(TVb{i}, 2))];
%                 end
%                 pclviewer(pcl_model);
                
                
                
                
                
                
                
                
                
                
                
                R0 = cell2mat(T(1,1,end))'*cell2mat(T(2,1,end));
                T0 = cell2mat(T(1,1,end))'*cell2mat(T(2,2,end)) - cell2mat(T(1, 2, end));

                
%                 R0 = cell2mat(T(1,1,end))'*cell2mat(T(2,1,end));
%                 T0 = cell2mat(T(1,1,end))'*cell2mat(T(2,2,end)) - cell2mat(T(1, 2, end));
%                 V2t = bsxfun(@plus, R0*V{2}, T0);
%                 pclviewer([[V{1} ; C{1}] [ V2t ; C{2}] ]);
%                 
%                 params.type = 2;
%                 params.Assigned = TAssigned;
%                 params.K = K;
%                 params.interval = 3;
%                 params.marker = marker;
%                 params.markerSize = markerSize;
%                 params.clrmap = clrmap;
%                 params.strIdx = strIdx;
%                 params.view = [];%[40 54];
%                 params.pause = 1;
%                 params.TXQ = TXQ;
%                 h = figure;
%                 [TV] = drawTransformation(V, T, params);


                T1 = R1 * T0 + T1;
                R1 = R1 * R0;
                

                curModel = netInfo(curModel);
            end
        end

        fprintf(sprintf('\n'));
    end
    
    trans.R = R1;
    trans.T = T1;
    
    disp(sprintf('put... %s', data_set.name{i}));
    transMatToOrigin(data_set.name{i}) = trans;
end



Rv = [];
Rn = [];
Rc = [];
RcSeg = [];

for i = 1:length(dataSetSampled)
    cData = idToData(i);%(data_set.name{i});
    transData = transMatToOrigin(data_set.name{i});
    
%     transData.R
%     transData.T

%     validIdx = find(cData.radius < 0.03);
%     Vt = bsxfun(@plus, transData.R*cData.vbf2(:,validIdx), transData.T);
%     Nt = bsxfun(@plus, transData.R*cData.nbf(:,validIdx), transData.T);
%     Ct = cData.cb(:,validIdx);

    Vt = bsxfun(@plus, transData.R*cData.vbf2, transData.T);
    Nt = bsxfun(@plus, transData.R*cData.nbf, transData.T);
    Ct = cData.cb;
%     Vt = cData.vs;
%     vaa = [vaa [ Vt ; cData.cs]];

    
    nameParts = strsplit(data_set.name{i}, '.');
    fNameOnly = nameParts{1};
    fileName = sprintf('%s', fNameOnly);
    %Save original size
    fcn_saveUniformSizeModel( Vt', Ct', Nt', dataOutRoot, fileName, 0 );

    Rv = [Rv ; Vt'];
    Rn = [Rn ; Nt'];
    Rc = [Rc ; Ct'];
    RcSeg = [RcSeg ; repmat(colors(i,:), size(Vt, 2), 1)];
end
pclviewer([Rv Rc]');
% pclviewer([Rv(1:10:end, :) Rc(1:10:end, :)]');
% pclviewer(vaa);
% 
fileName = 'all';
fcn_saveUniformSizeModel( Rv, Rc, Rn, dataOutRoot, fileName, 0 );
% 
fileName = 'segment';
fcn_saveUniformSizeModel( Rv, RcSeg, Rn, dataOutRoot, fileName, 0 );


params_desc.normalRadius=0.2;
params_desc.searchRadius=1.2;
params_desc.searchK = 0;
for i = 1:length(dataSetSampled)
    disp(sprintf('doing feature extraction... %d/%d', i, length(dataSetSampled)));
    cData = idToData(i);%(data_set.name{i});
    transData = transMatToOrigin(data_set.name{i});
    
%     transData.R
%     transData.T
    
    v1 = bsxfun(@plus, transData.R*cData.vsf2, transData.T)';
    n1 = bsxfun(@plus, transData.R*cData.nsg, transData.T)';
    c1 = cData.cs';

%     v1 = data.vs{i}';
%     c1 = data.cs{i}';
%     n1 = data.ns{i}';


    nameParts = strsplit(data_set.name{i}, '.');
    fNameOnly = nameParts{1};
%     fileName = sprintf('%s_ts', fNameOnly);
    %Save original size
%     fcn_saveUniformSizeModel( v1, c1, n1, dataRoot, fileName, 0 );

    
    [keypoints,~,~] = uniformSubSample(v1, 5, c1);
    descriptors = FPFHSDescriptor(v1, c1, keypoints, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);

    feature = [];
    feature.name = data_set.name{i};
    feature.keypoints = keypoints;
    feature.descriptors = descriptors;
    feature.v = v1;
    feature.c = c1;
    feature.n = n1;
    nameParts = strsplit(data_set.name{i}, '.');
    fNameOnly = nameParts{1};
    save(sprintf('%s/%s.mat', dataOutRoot, fNameOnly), 'feature');
end







% density = 50;
% %Save reduced size
% fileName = 'merged_ref_flat_ground3';
% fcn_saveUniformSizeModel( Rv, Rc, Rn, dataRoot, fileName, 50 );




return;









return;
idx = [5, 7];

V = {};
N = {};
C = {};

for i=idx
    V{end+1} = idToData(i).vs;
    C{end+1} = idToData(i).cs;
    N{end+1} = idToData(i).ns;
end

V = V';
C = C';
N = N';

% V = data.vs(idx);
% N = data.ns(idx);
% C = data.cs(idx);

[R,t,X,S,a,pk,T,TAssigned, TXQ] = joint_align(V,N,C);

params.type = 2;
params.Assigned = TAssigned;
params.K = size(TXQ{1,end}, 2);
params.interval = 5;
% params.marker = marker;
% params.markerSize = markerSize;
% params.clrmap = clrmap;
% params.strIdx = strIdx;
params.view = [];%[40 54];
params.pause = 1;
params.TXQ = TXQ;
h = figure;
[TV] = drawTransformation(V, T, params);


V = {};
N = {};
C = {};

for i=idx
    V{end+1} = idToData(i).vb;
    C{end+1} = idToData(i).cb;
    N{end+1} = idToData(i).nb;
end

V = V';
C = C';
N = N';

TV = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V,T(:,1,end),T(:,2,end),'uniformoutput',false);

vaa = [];
for i = 1:length(TV)
    vaa = [vaa [ TV{i} ;C{i}]];
end
% pclviewer([V{1} ;C{1}]);
pclviewer(vaa);








%%%%%%%%
[R2,t2,X2,S2,a2,pk2,T22,TAssigned2, TXQ2] = joint_align(V,N,C);

params2.type = 2;
params2.Assigned = TAssigned2;
params2.K = K;
params2.interval = 5;
params2.marker = marker;
params2.markerSize = markerSize;
params2.clrmap = clrmap;
params2.strIdx = strIdx;
params2.view = [];%[40 54];
params2.pause = 1;
params2.TXQ = TXQ2;
h = figure;
[TV2] = drawTransformation(V, T22, params2);



vaa = [];
for i = 1:length(TV2)
    vaa = [vaa [ TV2{i} ;C{i}]];
end
% pclviewer([V{1} ;C{1}]);
pclviewer(vaa);





































pclviewer([[V{1} ; C{1}] [ V{2} ; C{2}] ]);

R0 = cell2mat(T(1,1,end))'*cell2mat(T(2,1,end));
T0 = cell2mat(T(1,1,end))'*cell2mat(T(2,2,end)) - cell2mat(T(1, 2, end));
V2t = bsxfun(@plus, R0*V{2}, T0);
pclviewer([[V{1} ; C{1}] [ V2t ; C{2}] ]);


TV = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V,T(:,1,end),T(:,2,end),'uniformoutput',false);
pclviewer([[TV{1} ; C{1}] [ TV{2} ; C{2}] ]);

% vaa = [];
% for i = idx
%     vaa = [vaa [ TV{i} ;C{i}]];
% end
% % pclviewer([V{1} ;C{1}]);
% pclviewer(vaa);




