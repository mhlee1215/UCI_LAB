function [ vs, index, varargout ] = sampling( v, varargin )
%SAMPLING Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
% p.addOptional('ntimes',1);

p.addParameter('type','uniform');
p.addParameter('density', 10);
p.addParameter('binCapacity',0);
p.addParameter('color', []);
p.addParameter('normal', []);
p.addParameter('minNormalStd',0.05);
p.addParameter('minGroupSize',200);
p.parse(varargin{:});

varargout = {};

if strcmp(p.Results.type, 'uniform') 
    [vs, index, cs, ns] = uniformSubSample(v, p.Results.density, p.Results.color, p.Results.normal);
%     varargout{end+1} = sIndex;
    vs = vs';
    varargout{end+1} = cs';
    varargout{end+1} = ns';
elseif strcmp(p.Results.type, 'octree') 
    OT = OcTree(v, 'binCapacity', p.Results.binCapacity);
    
    
    PointBins = OT.PointBins;
    BinParents = OT.BinParents;
    BinDepths = OT.BinDepths;
%     
%     for i=1:max(OT.BinDepths)
%         leafBinMark(OT.BinParents(find(OT.BinDepths == i))) = 0;
%     end
%     leafBinIdx = find(leafBinMark);
%     
%     childPointIdx = find(ismember(OT.PointBins, leafBinIdx));
%     validPointBins = OT.PointBins(childPointIdx);
    minNormalStd = p.Results.minNormalStd;%0.05;
    minGroupSize = p.Results.minGroupSize;%100;
    %Leaf node statistics
    
    
    nextWorkMark = ones(1, OT.BinCount);%find(OT.BinDepths == max(OT.BinDepths));
%     eliminatedMark = zeros(1, OT.BinCount);%find(OT.BinDepths == max(OT.BinDepths));
%     for i=max(OT.BinDepths):-1:1
%         nextWorkMark(OT.BinParents(find(OT.BinDepths == i))) = 0;
%         
%     end
    
    
    for i=max(OT.BinDepths):-1:1
        
        %Remove parents which childeren are already taken
%         binIdxChildAlreadyTaken = find((leafBinMark == 0) .* (BinDepths == i));
%         if ~isempty(binIdxChildAlreadyTaken)
%             leafBinMark(BinParents(binIdxChildAlreadyTaken)) = 0;
% %             leafRemovalPointIdx = find(ismember(PointBins, binIdxChildAlreadyTaken));
% % %             PointBins(leafRemovalPointIdx) = PointBins(BinParents(leafRemovalPointIdx));
% %             BinParents(leafRemovalPointIdx) = BinParents(BinParents(leafRemovalPointIdx));
% %             BinDepths(leafRemovalPointIdx) = BinDepths(leafRemovalPointIdx)-1;
%         end
        %Get statics
        binIdxAlivePerDepth = find((nextWorkMark == 1) .* (BinDepths == i));
        if isempty(binIdxAlivePerDepth)
            continue;
        end
        leafIdx = find(ismember(PointBins, binIdxAlivePerDepth));
        outSize = [max(binIdxAlivePerDepth) 1];
        groupSize = accumarray(PointBins(leafIdx), ones(length(leafIdx), 1), outSize, @sum);
%         groupSizeContainer = zeros(max(binIdxAlivePerDepth), 1);
%         groupSizeContainer(1:length(groupSize)) = groupSize;
        %Eleminate node when it does not satisfy constraints         
        markOverConstraint = (groupSize(binIdxAlivePerDepth) < minGroupSize);
        if ~isempty(p.Results.normal)
            normal = p.Results.normal;
            sNXStd = accumarray(PointBins(leafIdx), normal(leafIdx,1), outSize, @std);
            sNYStd = accumarray(PointBins(leafIdx), normal(leafIdx,2), outSize, @std);
            sNZStd = accumarray(PointBins(leafIdx), normal(leafIdx,3), outSize, @std);
            sNStd = (sNXStd+sNYStd+sNZStd)./3;
            markOverConstraint = min(markOverConstraint + (sNStd(binIdxAlivePerDepth) < minNormalStd), 1);
        end
%         markOverConstraint = markOverConstraint(min(PointBins(leafIdx)):end);
        
        if sum(markOverConstraint)==0
            continue;
        end
        
        binIdxOverStd = binIdxAlivePerDepth(find(markOverConstraint));%+min(PointBins(leafIdx))-1;
        
        childSizeAll = accumarray((BinParents+1)', ones(length(BinParents), 1), [max(BinParents)+1 1], @sum);
        childSizeAll = childSizeAll(2:end);
        max(childSizeAll)
        
        childSize = accumarray((BinParents(binIdxOverStd)+1)', ones(length(binIdxOverStd), 1), [max(BinParents)+1 1], @sum);
        childSize = childSize(2:end);
%         childSize = childSize(min(BinParents(binIdxOverStd)):end);
%         binIdxAllOverParent = find(childSize>=7);
        binIdxAllOverParent = find((childSize > 0) .* ((childSize-childSizeAll) == 0) );
        
        childIdx = ismember(BinParents(binIdxOverStd), binIdxAllOverParent);
        binIdxValidOverStd = binIdxOverStd(childIdx);
        
%         childIdx = ismember(BinParents, binIdxAllOverParent);
%         binIdxValidOverStd = find(childIdx);
        
%         childSizeAll = accumarray((BinParents+1)', ones(length(BinParents), 1)', [], @sum);
        
        nextWorkMark(binIdxAllOverParent) = 1;
%         binIdxUnderStd = find(~markOverConstraint)+min(PointBins(leafIdx))-1;
%         nextWorkMark(BinParents(binIdxUnderStd)) = 0;
        
        
%         binIdxPerDepth = find((BinDepths == i));
        fromIdx = find(ismember(PointBins, binIdxValidOverStd));
%         toIdx = ismember(PointBins, BinParents(binIdxValidOverStd));
        PointBins(fromIdx) = BinParents(PointBins(fromIdx));
        
        vIdx = find(BinParents(binIdxAllOverParent) > 0);
        BinParents(binIdxAllOverParent(vIdx)) = BinParents(BinParents(binIdxAllOverParent(vIdx)));
        BinDepths(binIdxAllOverParent) = BinDepths(binIdxAllOverParent)-1;
    end    
    
    length(unique(OT.PointBins))
    length(unique(PointBins))
    
%      leafBinIdx = find(leafBinMark);
%      childPointIdx = find(ismember(PointBins, leafBinIdx));
    validPointBins = PointBins;%(childPointIdx);
    
    sXMean = accumarray(validPointBins, v(:,1), [], @mean);
    sYMean = accumarray(validPointBins, v(:,2), [], @mean);
    sZMean = accumarray(validPointBins, v(:,3), [], @mean);
    vs = [sXMean(:) sYMean(:) sZMean(:)];
    index = PointBins;
    
%     varargout{end+1} = BinDepths;
    
    
    if ~isempty(p.Results.color)
        color = p.Results.color;
        sRMean = accumarray(validPointBins, color(:,1), [], @mean);
        sGMean = accumarray(validPointBins, color(:,2), [], @mean);
        sBMean = accumarray(validPointBins, color(:,3), [], @mean);
        cs = [sRMean(:) sGMean(:) sBMean(:)];
        varargout{end+1} = cs;
    end
    if ~isempty(p.Results.normal)
        normal = p.Results.normal;
        sXMean = accumarray(validPointBins, normal(:,1), [], @mean);
        sYMean = accumarray(validPointBins, normal(:,2), [], @mean);
        sZMean = accumarray(validPointBins, normal(:,3), [], @mean);
        ns = [sXMean(:) sYMean(:) sZMean(:)];
        varargout{end+1} = ns;
    end
end



end

