function [ visible, visibleX, existX, visibleXBack ] = getCameraVisiblePoints2( V, pSet, params )
%GETVISIBLEPOINTS Summary of this function goes here
%   Detailed explanation goes here

X = [];
visibleX= [];
%In FOV but blocked by something
visibleXBack= [];
existX = [];
% srcId = 1;
% targetId = 1;

srcId = params.srcId;
% vertices = V;

% V = V;
TV = params.TV;

if nargin > 2
    emR = params.R;
    emT = params.T;
%     vertices = emR'*bsxfun(@minus, vertices, emT);
end


if isfield(params, 'X')
    X = params.X;
    visibleX = zeros(size(X, 2), 1);
    visibleXBack = zeros(size(X, 2), 1);
    existX = zeros(size(X, 2), 1);
end

camInterval = 50;
if isfield(params, 'camInterval')
    camInterval = params.camInterval;
end

colors = params.colors;
visParam = params.visParam;
visible = zeros(size(V, 2), 1);

% pSet = Pose_Set{srcId};
% figure; hold on;
%[visiblePtInds, inVisiblePtInds, inVisibleVector]=HPR(p,C,param)

for pi = 1:camInterval:length(pSet)
    sprintf('%d/%d\n', pi, length(pSet))
    t = pSet(pi).t;
    r = pSet(pi).r;

    
    c_pos = t;
    c_pos_em = emR*c_pos + emT;
    
%     v2 = bsxfun(@minus, r'*vertices, t);
    v2 = r'*bsxfun(@minus, V, t);
    front = v2(3,:)>0;
    v2p = v2 ./ repmat(v2(3,:), 3, 1);
    curVis = front .* (v2p(1,:) < 0.55) .* (v2p(1,:) > -0.55) ...
                        .* (v2p(2,:) < 0.46) .* (v2p(2,:) > -0.46);
    
   
    
%     curVisIdx = find(curVis);
%     newVisIdx = curVisIdx(~inVisibleVector(1:end-lenX));
%     size(newVisIdx)
%     curVis = curVis .* 0;
%     curVis(newVisIdx) = 1;
    
    visible = max(visible, curVis');
   
    curVisX = [];
    if ~isempty(X)
        X2 = r'*bsxfun(@minus, X, t);
        frontX = X2(3,:)>0;
        X2 = X2 ./ repmat(X2(3,:), 3, 1);
        curVisX = frontX .* (X2(1,:) < 0.55) .* (X2(1,:) > -0.55) .* (X2(2,:) < 0.46) .* (X2(2,:) > -0.46);
        
        
        
        curVisColor = colors(:, find(curVis));
        %For each emRvisible X (within FOV)
        curVisXPoint = X(:, find(curVisX));
        visibleX = max(visibleX, curVisX');
        pinXVisAll = {};
%         figure; 
        for xIdx = find(curVisX)
            if existX(xIdx) == 1
                continue;
            end
%             find((abs(X(1,:) - -2.115) < 0.01) .* (abs(X(2,:) - -1.116) < 0.01))
%             pIdxinX = [];
%             for i=1:length(params.clustAssign)
%                 pIdxinX = find(params.clustAssign==xIdx);
%             end
            
            
            
%             figure;
            pinX = [];%vertices(:,pIdxinX);
%             
            %Collect every points which belongs to cluster xIdx
            for i=1:length(params.clustAssign)
%                 for xIdx = find(curVisX)
                    
                    pinXSub = params.TV{i}(:, find(params.clustAssign{i}==xIdx));
                    pinX = [pinX pinXSub];
%                 end
            end
            
            %points within a cluster. Should be removed from compPoints
            %sets. Otherwise, it can be blocked by themselves.
            pinXSelfVec = params.clustAssign{srcId}==xIdx;
            
            
            %Get compPoints itself. 
%             clusterPoints = TV{srcId}(:, find(pinXSelfVec));
            
            %If a cluster has lots of points which belongs to it, it is
            %exist
            if sum(pinXSelfVec) > 0
                existX(xIdx) = 1;
                continue;
            end
            
%             Get compPoints from visible point set without cluster itself
            compPoints = TV{srcId}(:, find(curVis.*~pinXSelfVec));
            
            if isempty(pinX)
                continue;
            end
            
            maxCheckP = 15;
            randIdx = randperm(size(pinX, 2));
            pinX = pinX(:, randIdx(1:min(maxCheckP, size(pinX, 2))));
            pinXhasBack = zeros(size(pinX, 2), 1);
            for pId = 1:size(pinX, 2)
                pX = pinX(:, pId);
%                 pX = [-2.12 -1.066 -0.4755]';
                %get Eucledian distance from camera to point for check
                distFromCam = norm(pX - c_pos_em);
                compLine = [c_pos_em' pX'-c_pos_em'];
                
                %Get distance from line to points
                d = distancePointLine3d(compPoints', compLine);
                %Get projected points from points to line
                projP = projPointOnLine3d(compPoints', compLine);
                %Get distance from camera center to projected points
                distToPoints = sqrt(sum(bsxfun(@minus, projP', c_pos_em).^2));
                %Count how many points are behind of the point (pX)
                inPointVec = (distFromCam + 0.7 < distToPoints) .* ( d' < 0.03) ;
                inPointNum = sum(inPointVec);
                inPointIdx = find(inPointVec);
                
                if inPointNum > 0
                    pinXhasBack(pId) = 1;% = [ pinXVis 1];
                else
                    pinXhasBack(pId) = 0;% = [ pinXVis 0];
                end
%                   pclviewer([compPoints c_pos_em; ones(size(projP')) [1 0 0]']);
%                  pclviewer([compPoints projP' c_pos_em ; [(d./max(d))' ; zeros(2, length(d))] ones(size(projP')) [0 1 0]']);
%                 pclviewer([compPoints(:, inPointIdx) projP' c_pos_em ]);
%                 pclviewer([compPoints(:, :) projP' c_pos_em pX ;  [(d./max(d))' ; zeros(2, length(d))] ones(size(projP')) [0 1 0]' [0 0 1]']);
%                 pclviewer([compPoints(:, :) c_pos_em pX curVisXPoint;  [(d./max(d))' ; zeros(2, length(d))] [0 1 0]' [0 0 1]' ones(size(curVisXPoint))]);
%                 
%                 pnts = [compPoints(:, inPointIdx) c_pos_em pX curVisXPoint];
%                 clrs = [[(d(inPointIdx)./max(d))' ; zeros(2, length(d(inPointIdx)))] [0 1 0]' [0 0 1]' repmat([0 0.5 1]', 1, size(curVisXPoint, 2))];
                
                pnts = [compPoints(:, :) c_pos_em pX curVisXPoint];
                clrs = [[(d(:)./max(d))' ; zeros(2, length(d(:)))] [0 1 0]' [0 0 1]' repmat([0 0.5 1]', 1, size(curVisXPoint, 2))];
                
%                 figure; clf;
%                 scatter3(pnts(1,:), pnts(2,:), pnts(3,:), 8, clrs(:, :)', 'filled'); hold on;
%                 scatter3(pX(1,:), pX(2,:), pX(3,:), 35, 'g', 'filled'); 
%                 scatter3(X(1,xIdx), X(2,xIdx), X(3,xIdx), 35, 'm', 'filled'); 
%                 scatter3(pinX(1,find(pinXhasBack)), pinX(2,find(pinXhasBack)), pinX(3,find(pinXhasBack)), 35, 'b', 'filled'); 
%                 scatter3(pinX(1,find(~pinXhasBack)), pinX(2,find(~pinXhasBack)), pinX(3,find(~pinXhasBack)), 35, 'r', 'filled'); 
%                 
%                 
%                 line([c_pos_em(1) pX(1)], [c_pos_em(2) pX(2)], [c_pos_em(3) pX(3)], 'linewidth', 2, 'color', 'm'); 
%                 title(sprintf('xIdx :%d', xIdx));
%                 axis equal;
%                 view(-206, 22);
%                 pause(0.1);
            end
            
            pinXVisAll{end+1} =  pinXhasBack;
            
            
            if mean(pinXhasBack) > 0.7
                existX(xIdx) = 0;
            else
                visibleXBack(xIdx) = 1;
            end
            
%                 pclviewer([X(:, find(curVisX.*existX')) X(:, find(curVisX.*~existX')) ; ...
%                     repmat([1 0 0]', 1, sum(curVisX.*existX')) repmat([0 1 0]', 1, sum(curVisX.*~existX'))]);
            
        end
    end
    
    
    sprintf('hi\n');
    
    
    if ~isempty(X)
        [~,~,inVisibleVector]=HPR([TV{srcId}(:, find(curVis)) X(:, find(curVisX))]', c_pos_em',visParam);
        lenX = sum(curVisX);
        curVisIdxX = find(curVisX);
        newVisIdxX = curVisIdxX(~inVisibleVector(end-lenX+1:end));
%         size(newVisIdxX)
        curVisX = curVisX .* 0;
        curVisX(newVisIdxX) = 1;
        visibleX = max(visibleX, curVisX');
    end

    %%Visualization
    
%     corners = [ 0.55 0.44 1 ; -0.55 0.44 1 ;-0.55 -0.44 1; 0.55 -0.44 1]';
%     corners = corners .*3;
%     c_world2 = bsxfun(@plus, r*corners, t);
%     clf;
%     
% %     scatter3(vertices(1, find(curVis)), vertices(2, find(curVis)), vertices(3, find(curVis)), 8, colors(:, find(curVis))', 'filled');
%     
%     colors2 = colors;
%     colors2(:, find(curVis)) = bsxfun(@plus, colors2(:, find(curVis)), [0.5 0 0]');
%     
%     scatter3(vertices(1, :), vertices(2, :), vertices(3, :), 8, colors2(:, :)', 'filled');
%     
%     hold on;
%     scatter3(c_pos(1), c_pos(2), c_pos(3), 10, [1 0 0], 'filled');
%     
%     scatter3(0, 0, 0, 15, [1 0 0], 'filled');
%     
% %     for li=1:4
% %         line([c_pos(1) c_world(1, li)]', [c_pos(2) c_world(2, li)]', [c_pos(3) c_world(3, li)]', 'linewidth', 2);    
% %         line([c_world(1, li) c_world(1, mod(li, 4)+1)]', [c_world(2, li) c_world(2, mod(li, 4)+1)]', [c_world(3, li) c_world(3, mod(li, 4)+1)]', 'linewidth', 2, 'color', 'g'); 
% %     end
%     
%     for li=1:4
%         line([c_pos(1) c_world2(1, li)]', [c_pos(2) c_world2(2, li)]', [c_pos(3) c_world2(3, li)]', 'linewidth', 2, 'color', 'k');    
%         line([c_world2(1, li) c_world2(1, mod(li, 4)+1)]', [c_world2(2, li) c_world2(2, mod(li, 4)+1)]', [c_world2(3, li) c_world2(3, mod(li, 4)+1)]', 'linewidth', 2, 'color', 'm'); 
%     end
%     
% %     li=1;
% %     line([c_world(1, li) c_world(1, mod(li, 4)+1)]', [c_world(2, li) c_world(2, mod(li, 4)+1)]', [c_world(3, li) c_world(3, mod(li, 4)+1)]', 'linewidth', 1, 'color', 'g');    
%     
%     view(182, -53);
% %     view(-180, -53);
% %     view(63, -55);
%     axis equal;
% %     view(-156, 34);
% %     view(186, -11);
%     axis([-5 5 -3 3 -5 5]);
%     pause(0.1);
    
end


end

