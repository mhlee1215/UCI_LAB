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

% clc
% close all
% clear all
% % g = gpuDevice(1);
% % reset(g);
% 
% addpath(genpath('libs'));
% % addpath(genpath('libs2'));
% % run('libs2/vlfeat-0.9.20-bin/vlfeat-0.9.20/toolbox/vl_setup');
% 
% 
% % data3DFilePath = '/home/mhlee/Kinect_Logs/2016-02-03.00.klg.ply';
% % [tri, pts, data, comments] = ply_read(data3DFilePath, 'tri');
% % 
% % v = [data.vertex.x data.vertex.y data.vertex.z];
% % 
% % 
% % figure; fv_alligned.Faces = tri;
% % fv_alligned.Vertices = v;
% % patch(fv_alligned,'FaceColor',       [0.8 0.8 1.0], ...
% %          'EdgeColor',       'none',        ...
% %          'FaceLighting',    'gouraud',     ...
% %          'AmbientStrength', 0.15);
% % camlight('headlight');
% % material('dull');
% % axis equal;
% % hold on;
% 
% %data_set = [1:6 8:14];%[1 2 3 4];
% dataRoot = '/home/mhlee/data_from_odroid/complete';
% % data_set = {'2016-02-09.01','2016-02-09.02','2016-02-09.03'};
% % data_set = {'2016-04-29.03', '2016-04-29.05', '2016-04-29.06', '2016-04-29.00','2016-04-29.08'};
% %Wrong init pos : '2016-04-29.06'
% % data_set = {'2016-04-29.00', '2016-04-29.01', '2016-04-29.02', '2016-04-29.03', '2016-04-29.05', '2016-04-29.08'};%, '2016-04-29.05'};%, '2016-04-29.00','2016-04-29.08'};
% data_set = {'MAC_TEST-2016-07-17_23_11'};
% d_size = length(data_set);
% 
% 
% [dataSet, poseSet] = loadEFDataset( dataRoot, data_set);
% dataSetSampled = loadUniformSampling(dataSet, 10);
% 
% V = {};
% N = {};
% C = {};
% for i=1:length(dataSetSampled)
%     V{i} = dataSetSampled{i}.v;
%     N{i} = dataSetSampled{i}.n;
%     C{i} = dataSetSampled{i}.c;
% end
% V = V';
% N = N';
% C = C';

function [] = Annotation3DGui(params)

Vertex = params.Vertex;
Color = params.Color;
SegVertex = params.SegVertex;
SegColor = params.SegColor;
Seg = params.Seg;

global annotationResult;
annotationResult = {};

global prevSelectedObjectIndex;
prevSelectedObjectIndex = -1;

global selectedObjectIndex;
selectedObjectIndex = 1;

global isFigTop;
isFigTop = true;
global prevGroup;
prevGroup = 0;

global visibleGroupLeftIdx;
visibleGroupLeftIdx = ones(length(Seg.group), 1);

global visibleFigIdx;
visibleFigIdx = ones(length(Vertex), 1);

global visibleFigIdx2;
visibleFigIdx2 = {};%ones(length(Vertex), 1);

global seg;
seg = Seg;

global V;
global C;


global S;

% [ segmentedPoints, coloredCloud ] = Segmentation(Vertex{i}', Color{i}', 1500, '');

V = Vertex;%coloredCloud(1:3, :);
C = Color;%coloredCloud(4:6, :);
SV = SegVertex;
SC = SegColor;

bWidth = 130;

paddingFig = 50;
height = 580;
% ws = [1280 680];
fSize = [height-paddingFig*2 height-paddingFig*2];
width = 1280;
ws = [width height];

SCR = get(0,'Screensize');  % Get screensize.
S.fh = figure('numbertitle','off',...
              'menubar','none',...
              'units','pixels',...
              'position',[SCR(3)/2-ws(1)/2 ,SCR(4)/2-ws(2)/2 , ws(1), ws(2)],...
              'name','Point Cloud Labeling Tool',...
              'resize','off');
S.ax = axes('units','pixels',...
            'position',[paddingFig paddingFig fSize(1) fSize(2)]);
hold on;    
        
S.axt2 = {};
S.axt_seg2 = {};

S.ax2 = axes('units','pixels',...
            'position',[paddingFig*3+fSize(1)+bWidth paddingFig fSize(1) fSize(2)]);
hold on;    
% 
% axis equal;
grid(S.ax);
grid(S.ax2);
axis(S.ax, 'equal');
axis(S.ax2, 'equal');

for i=1:length(Seg.group)
    S.axt{i} = scatter3(S.ax, V(1,Seg.group{i})', V(2,Seg.group{i})', ...
        V(3,Seg.group{i})', 8, C(:,Seg.group{i})', 'filled'); 
    S.axt_seg{i} = scatter3(S.ax, SV(1,Seg.group{i})', SV(2,Seg.group{i})', ...
        SV(3,Seg.group{i})', 8, SC(:, Seg.group{i})', 'filled');    
    set(S.axt_seg{i}, 'Visible', 'off');
    
    S.axt2{i} = scatter3(S.ax2, V(1,Seg.group{i})', V(2,Seg.group{i})', ...
        V(3,Seg.group{i})', 8, C(:,Seg.group{i})', 'filled'); 
    set(S.axt2{i}, 'Visible', 'off');
    S.axt_seg2{i} = scatter3(S.ax2, SV(1,Seg.group{i})', SV(2,Seg.group{i})', ...
        SV(3,Seg.group{i})', 8, SC(:, Seg.group{i})', 'filled');    
    set(S.axt_seg2{i}, 'Visible', 'off');
    
    visibleFigIdx2{i} = zeros(length(Vertex), 1);
end



% axis('manual');
% range(reshape([min(Vertex') ;max(Vertex')], 1, 6));





paddingTop = 15;
padding = 10;
bPos = [fSize(1)+paddingFig*2 ws(2)];

bPos = bPos - [0 paddingTop];          


bSize = [bWidth 40];
S.pb1 = uicontrol('style','pushbutton',...
                  'units','pixels',...
                  'position',[bPos(1) bPos(2)-bSize(2) bSize(1) bSize(2)],...
                  'string','Flip',...
                  'fontsize',12, ...
                  'callback',{@flip_call});

bPos = bPos - [0 bSize(2)+padding];              
bSize = [bWidth 180];


defaultObjects = {'floor';'wall';'chair1';'chair2';'desk1';'desk2'};
for i=1:length(defaultObjects)
    annotationResult{i}.name = defaultObjects{i};
    annotationResult{i}.group = zeros(1, length(Seg.group));
end
S.ls = uicontrol('style','list',...
                 'unit','pix',...
                 'position',[bPos(1) bPos(2)-bSize(2) bSize(1) bSize(2)],...
                 'min',0,'max',2,...
                 'fontsize',12,...
                 'string', defaultObjects...
                 );         
bPos = bPos - [0 bSize(2)+padding];                           
bSize = [bWidth 30];             
S.ed = uicontrol('style','edit',...
                 'unit','pix',...
                 'position',[bPos(1) bPos(2)-bSize(2) bSize(1) bSize(2)],...
                 'fontsize',12,...
                 'string','New Object');
bPos = bPos - [0 bSize(2)+padding];                           
bSize = [bWidth 40];             
S.pb = uicontrol('style','push',...
                 'units','pix',...
                 'position',[bPos(1) bPos(2)-bSize(2) bSize(1) bSize(2)],...
                 'fontsize',12,...
                 'string','Add Object',...
                 'callback',{@ed_call, Seg});
             
axis equal;
% S.ax2 = axes('units','pixels',...
%             'position',[50 50 ws(2)-100 ws(2)-100]);
% S.axt2 = scatter3(V(1,:)', V(2,:)', V(3,:)', 8, C', 'filled');
% 
% set(S.ax2, 'Visible', 'off');


cameratoolbar('Show');
set( gcf, 'menubar', 'figure' )
set(gcf,'toolbar','figure');
hold on;


S.updateLeftFig = @updateLeftFig;
S.updateRightFig = @updateRightFig;


set(S.ax, 'ButtonDownFcn', {@myCallbackClickA3DPoint, SV, C, Seg}); 
for i=1:length(Seg.group)
    set(S.axt{i}, 'ButtonDownFcn', {@myCallbackClickA3DPoint, SV, C, Seg}); 
    set(S.axt_seg{i}, 'ButtonDownFcn', {@myCallbackClickA3DPoint, SV, C, Seg}); 
end

set(S.ax2, 'ButtonDownFcn', {@myCallbackClickA3DPointRight, SV, C, Seg}); 
for i=1:length(Seg.group)
    set(S.axt2{i}, 'ButtonDownFcn', {@myCallbackClickA3DPointRight, SV, C, Seg}); 
    set(S.axt_seg2{i}, 'ButtonDownFcn', {@myCallbackClickA3DPointRight, SV, C, Seg}); 
end

set(S.ls, 'callback',{@list_call});

S.testFun = @testFun;       
end

function [] = ed_call(varargin)

global annotationResult;

% Callback for pushbutton, adds new string from edit box.
% S = varargin{3};  % Get the structure.
global S;
global seg;
Seg = seg;
% Seg = varargin{4};

oldstr = get(S.ls,'string'); % The string as it is now.
addstr = {get(S.ed,'string')}; % The string to add to the stack.
% The order of the args to cat puts the new string either on top or bottom.
set(S.ls,'str',{addstr{:},oldstr{:}});  % Put the new string on top -OR-
% set(S.ls,'str',{oldstr{:},addstr{:}});  % Put the new string on bottom.

annotationResult{end+1}.name = addstr;
annotationResult{end}.group = zeros(1, length(Seg.group));



end

function [] = flip_call(varargin)
global annotationResult;
global visibleGroupLeftIdx;
global selectedObjectIndex;
global S;
% S = varargin{3};  % Get the structure.

global isFigTop;
    if isFigTop == true
        %uistack(S.ax2, 'up', 1);
        for i=1:length(visibleGroupLeftIdx)
            if visibleGroupLeftIdx(i) == 1
                set(S.axt{i}, 'Visible', 'off');
                set(S.axt_seg{i}, 'Visible', 'on');
            end
        end
        
        visibleSegGroup = annotationResult{selectedObjectIndex}.group;
        for i = 1:length(visibleSegGroup)
            if visibleSegGroup(i) == 1
                set(S.axt2{i}, 'Visible', 'off');
                set(S.axt_seg2{i}, 'Visible', 'on');
            end
        end
        
%         set(S.axt, 'Visible', 'off');
%         set(S.ax2, 'Visible', 'On');
%         set(S.axt2, 'Visible', 'On');
    else
        %uistack(S.ax, 'up', 1);
         for i=1:length(visibleGroupLeftIdx)
            if visibleGroupLeftIdx(i) == 1
                set(S.axt{i}, 'Visible', 'on');
                set(S.axt_seg{i}, 'Visible', 'off');
            end
        end
        
        visibleSegGroup = annotationResult{selectedObjectIndex}.group;
        for i = 1:length(visibleSegGroup)
            if visibleSegGroup(i) == 1
                set(S.axt2{i}, 'Visible', 'on');
                set(S.axt_seg2{i}, 'Visible', 'off');
            end
        end
        
    end
    
    isFigTop = ~isFigTop;
end

function [] = list_call(varargin)
% Callback for pushbutton, adds new string from edit box.
global S;
% S = varargin{3};  % Get the structure.

global prevSelectedObjectIndex;
global selectedObjectIndex;

if prevSelectedObjectIndex ~= selectedObjectIndex
    prevSelectedObjectIndex = selectedObjectIndex;
else
    prevSelectedObjectIndex = -1;
end
selectedObjectIndex = get(S.ls, 'value');

S.updateRightFig();

% oldstr = get(S.ls,'string'); % The string as it is now.
% addstr = {get(S.ed,'string')}; % The string to add to the stack.
% % The order of the args to cat puts the new string either on top or bottom.
% set(S.ls,'str',{addstr{:},oldstr{:}});  % Put the new string on top -OR-
% set(S.ls,'str',{oldstr{:},addstr{:}});  % Put the new string on bottom.
end

function[] = updateLeftFig()

global visibleGroupLeftIdx;
global prevSelectedObjectIndex;
global selectedObjectIndex;
global seg;
global V;
global C;
global S;
global isFigTop;
disp(sprintf('updateLeftFig Called!, prev:%d, cur:%d', prevSelectedObjectIndex, selectedObjectIndex));

    for i=1:length(visibleGroupLeftIdx)
            if visibleGroupLeftIdx(i) == 1
                if isFigTop == false
                    set(S.axt{i}, 'Visible', 'off');
                    set(S.axt_seg{i}, 'Visible', 'on');
                else
                    set(S.axt{i}, 'Visible', 'on');
                    set(S.axt_seg{i}, 'Visible', 'off');
                end
            end
    end
        

end

function[] = updateRightFig()
global annotationResult;
global prevSelectedObjectIndex;
global selectedObjectIndex;
global seg;
global V;
global C;
global S;
global isFigTop;
disp(sprintf('updateRightFig Called!, prev:%d, cur:%d', prevSelectedObjectIndex, selectedObjectIndex));

    if prevSelectedObjectIndex ~= -1
        prevVisibleSegGroup = annotationResult{prevSelectedObjectIndex}.group;
        for i = 1:length(prevVisibleSegGroup)
            if prevVisibleSegGroup(i) == 1
                set(S.axt2{i}, 'Visible', 'off');
                set(S.axt_seg2{i}, 'Visible', 'off');

            end
        end
    end
    
    visibleSegGroup = annotationResult{selectedObjectIndex}.group;
    for i = 1:length(visibleSegGroup)
        if visibleSegGroup(i) == 1
            if isFigTop == false
                set(S.axt2{i}, 'Visible', 'off');
                set(S.axt_seg2{i}, 'Visible', 'on');
            else
                set(S.axt2{i}, 'Visible', 'on');
                set(S.axt_seg2{i}, 'Visible', 'off');
            end
        end
    end
end
% pointCloud = V{1};
% h = gcf;
% % plot3(pointCloud(1,:), pointCloud(2,:), pointCloud(3,:), 'c.'); 
% scatter3(V{i}(1,:)', V{i}(2,:)', V{i}(3,:)', 8, C{i}', 'filled');
% cameratoolbar('Show'); % show the camera toolbar
% hold on; % so we can highlight clicked points without clearing the figure
% 
% % set the callback, pass pointCloud to the callback function
% set(h, 'WindowButtonDownFcn', {@callbackClickA3DPoint, pointCloud}); 
% S.pb2 = uicontrol('style','push',...
%                  'units','pix',...
%                  'position',[50 50 180 80],...
%                  'fontsize',14,...
%                  'string','Hi2!');   
% 
% clickA3DPoint(V{1});


