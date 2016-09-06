function [] = Annotation3DGui(params)
init(params);
end

function [] = init(params)


density = 20;



global dataIdx;
dataIdx = 1;

global plyData;
plyData = params.plyData(density);
global segData;
segData = params.Seg(density);

global dataRange;
minRange = [9999 9999 9999];
for i=1:length(plyData.v)
    minRange = min(minRange, min(plyData.v{i}'));
end

maxRange = [-9999 -9999 -9999];
for i=1:length(plyData.v)
    maxRange = max(maxRange, max(plyData.v{i}'));
end

dataRange.min = minRange;
dataRange.max = maxRange;


global dataSize;
dataSize = length(plyData.v);

global annotationResult;
for i=1:dataSize
    annotationResult{i} = {};
end

global prevSelectedObjectIndex;
prevSelectedObjectIndex = -1;

global selectedObjectIndex;
selectedObjectIndex = 1;

global isFigTop;
isFigTop = true;
global prevGroup;
prevGroup = 0;

global S;

global visibleGroupLeftIdx;

visibleGroupLeftIdx = {};
for i=1:dataSize
    visibleGroupLeftIdx{i} = ones(length(segData{i}.group), 1);
end
%
global visibleFigIdx;
visibleFigIdx = {};
for i=1:dataSize
    visibleFigIdx{i} = ones(length(plyData.v{i}), 1);
    %     for j=1:length(segData{i}.group)
    %         visibleFigIdx{i}{j} = ones(length(length(segData{i}.group{j})), 1);
    %     end
    
end
%
global visibleFigIdx2;
visibleFigIdx2 = {};
for i=1:dataSize
    for j=1:length(segData{i}.group)
        visibleFigIdx2{i}{j} = zeros(length(segData{i}.group{j}), 1);
    end
end

% global seg;
% seg = Seg;

% global V;
% global C;




% d = params.plyData(density);
% s = params.Seg(density);
% s = s{1};
%
% Vertex = d.v;
% Color = d.c;
% SegVertex = s.SegVertex;
% SegColor = s.SegColor;
% Seg = s;%params.Seg;
date_set = params.date_set;
name_set = params.name_set;


% [ segmentedPoints, coloredCloud ] = Segmentation(Vertex{i}', Color{i}', 1500, '');

V = plyData.v{dataIdx};%coloredCloud(1:3, :);
C = plyData.c{dataIdx};%coloredCloud(4:6, :);
SV = segData{dataIdx}.SegVertex;
SC = segData{dataIdx}.SegColor;
Seg = segData{dataIdx};

bWidth = 240;
topPanelHeight = 250;
bottomPaneHeight = 100;
paddingFig = 50;
height = 900;
width = 1400;
% ws = [1280 680];
fWidth = height-bottomPaneHeight-paddingFig*2-topPanelHeight;
fSize = [fWidth fWidth];

ws = [width height];

SCR = get(0,'Screensize');  % Get screensize.
S.fh = figure('numbertitle','on',...
    'menubar','none',...
    'units','pixels',...
    'position',[SCR(3)/2-ws(1)/2 ,SCR(4)/2-ws(2)/2 , ws(1), ws(2)],...
    'name','Point Cloud Labeling Tool',...
    'resize','off');
S.ax = axes('units','pixels',...
    'position',[paddingFig paddingFig+bottomPaneHeight fSize(1) fSize(2)]);
% hold on;



% TAssigned = params.TAssigned;
% unObservedAllX  = params.unObservedAllX;
% SEG_color =  params.SEG_color;
% TV = params.TV;

padding = 10;
bPos = [paddingFig paddingFig*1.5+bottomPaneHeight+fSize(2)+topPanelHeight];
bSize = [bWidth-50 150];
defaultObjects = {'floor';'wall';'chair1';'chair2';'desk1';'desk2'};
for i=1:length(defaultObjects)
    annotationResult{dataIdx}{i}.name = defaultObjects{i};
    annotationResult{dataIdx}{i}.group = zeros(1, length(Seg.group));
end

S.ls = uicontrol('style','list',...
    'unit','pix',...
    'position',[bPos(1) bPos(2)-bSize(2) bSize(1) bSize(2)],...
    'min',0,'max',2,...
    'fontsize',12,...
    'string', defaultObjects...
    );
bPos = bPos - [0 bSize(2)+padding];
bSize = [bWidth-50 30];
S.ed = uicontrol('style','edit',...
    'unit','pix',...
    'position',[bPos(1) bPos(2)-bSize(2) bSize(1) bSize(2)],...
    'fontsize',12,...
    'string','New Object');
bPos = bPos - [0 bSize(2)+padding];
bSize = [bWidth-50 30];
S.pb = uicontrol('style','push',...
    'units','pix',...
    'position',[bPos(1) bPos(2)-bSize(2) bSize(1) bSize(2)],...
    'fontsize',12,...
    'string','Add Object',...
    'callback',{@ed_call, Seg});

S.axSubt = {};
S.axSubt_seg = {};
for ii=1:4
    S.axSubt{ii} = {};
    S.axSubt_seg{ii} = {};
    S.axSub(ii) = axes('units','pixels',...
        'position',[bWidth+50+paddingFig+(ii-1)*(topPanelHeight-paddingFig/2+paddingFig/2) paddingFig*2+bottomPaneHeight+fSize(2) topPanelHeight-paddingFig/2 topPanelHeight-paddingFig/2]);
    
    %     hold on;
    %     for src = ii
    % %         h=figure; hold on;
    %         for i=src%1:length(TV)
    %             clustAssin = cell2mat(TAssigned(i,end));
    % %             observedVec = unObservedAllX(src, clustAssin);%max(observedAllX(src,clustAssin), observedAllX2(src,clustAssin));
    %             scatter3(S.axSub(ii), TV{i}(1,:)', TV{i}(2,:)', TV{i}(3,:)', 8, SEG_color(observedVec+1, :), 'filled');
    %         end
    %
    %         view(-184, -27);
    %
    %     end
    grid(S.axSub(ii));
    axis(S.axSub(ii), 'equal');
end



S.axt2 = {};
S.axt_seg2 = {};

S.ax2 = axes('units','pixels',...
    'position',[paddingFig*3+fSize(1)+bWidth paddingFig+bottomPaneHeight fSize(1) fSize(2)]);
% hold on;
%
% axis equal;

initFigure();
updateFigure();

% for i=1:length(Seg.group)
%     S.axt{i} = scatter3(S.ax, V(1,Seg.group{i})', V(2,Seg.group{i})', ...
%         V(3,Seg.group{i})', 8, C(:,Seg.group{i})', 'filled');
%     S.axt_seg{i} = scatter3(S.ax, SV(1,Seg.group{i})', SV(2,Seg.group{i})', ...
%         SV(3,Seg.group{i})', 8, SC(:, Seg.group{i})', 'filled');
%     set(S.axt_seg{i}, 'Visible', 'off');
%
%     S.axt2{i} = scatter3(S.ax2, V(1,Seg.group{i})', V(2,Seg.group{i})', ...
%         V(3,Seg.group{i})', 8, C(:,Seg.group{i})', 'filled');
%     set(S.axt2{i}, 'Visible', 'off');
%     S.axt_seg2{i} = scatter3(S.ax2, SV(1,Seg.group{i})', SV(2,Seg.group{i})', ...
%         SV(3,Seg.group{i})', 8, SC(:, Seg.group{i})', 'filled');
%     set(S.axt_seg2{i}, 'Visible', 'off');
%
% %     visibleFigIdx2{i} = zeros(length(Vertex), 1);
% end



% axis('manual');
% range(reshape([min(Vertex') ;max(Vertex')], 1, 6));





paddingTop = 15;
padding = 10;
bPos = [fSize(1)+paddingFig*2 ws(2)];

bPos = bPos - [0 paddingTop+topPanelHeight+paddingFig/2];


% bSize = [bWidth 40];
% S.pb1 = uicontrol('style','pushbutton',...
%                   'units','pixels',...
%                   'position',[bPos(1) bPos(2)-bSize(2) bSize(1) bSize(2)],...
%                   'string','RGB',...
%                   'fontsize',12, ...
%                   'callback',{@flip_call});
% bPos = bPos - [0 bSize(2)+padding];

bSize = [bWidth 20];
S.pp = uicontrol('style','pop',...
    'unit','pix',...
    'position',[bPos(1) bPos(2)-bSize(2) bSize(1) bSize(2)],...
    'backgroundc',get(S.fh,'color'),...
    'string',date_set,'value',1);
bPos = bPos - [0 bSize(2)+padding];

bSize = [bWidth 30];
S.bg = uibuttongroup('units','pix',...
    'pos',[bPos(1) bPos(2)-bSize(2) bSize(1) bSize(2)]);
% S.rd(1) = uicontrol(S.bg,...
%                     'style','rad',...
%                     'unit','pix',...
%                     'position',[10 7 110 30],...
%                     'string','Dynamic/Static',...
%                     'callback',{@radio_call,S});
S.rd(1) = uicontrol(S.bg,...
    'style','rad',...
    'unit','pix',...
    'position',[10 0 110 30],...
    'string','RGB',...
    'callback',{@radio_call,S});
S.rd(2) = uicontrol(S.bg,...
    'style','rad',...
    'unit','pix',...
    'position',[120 0 110 30],...
    'string','Segments',...
    'callback',{@radio_call,S});
bPos = bPos - [0 bSize(2)+padding];
% set(S.pb,'callback',{@pb_call,S}); % Set the callback, pass hands.

bSize = [bWidth 350];
% S.fh
% date = {};
% date{1} = '2016-02-02';
% date{2} = '2016-02-02';

v = {};%
for i=1:length(date_set)
    v{end+1} = true;
end

tData = [date_set' v'];
%
% load patients LastName Age Weight Height SelfAssessedHealthStatus
% PatientData = [LastName num2cell([Age Weight Height]) SelfAssessedHealthStatus];

t = uitable('Position', [bPos(1) bPos(2)-bSize(2) bSize(1) bSize(2)], 'Data', tData);
t.ColumnEditable = [false true];
% t.ColumnFormat = {[] [] [] [] {'Excellent', 'Fair', 'Good', 'Poor'}};
t.ColumnName = {'Date', 'visible?'};
t.ColumnWidth = {[130], 'auto'};
% t.ColumnFormat = {'string', 'logical'};
bPos = bPos - [0 bSize(2)+padding];



bPos = [paddingFig paddingFig];
bSize = [100 30];
S.pb = uicontrol('style','push',...
    'units','pix',...
    'position',[bPos(1) bPos(2)-bSize(2) bSize(1) bSize(2)],...
    'fontsize',12,...
    'string','< Prev',...
    'callback',{@go_prev, Seg});

bSize = [100 30];
bPos = [width-paddingFig-bSize(1) paddingFig];

S.pb = uicontrol('style','push',...
    'units','pix',...
    'position',[bPos(1) bPos(2)-bSize(2) bSize(1) bSize(2)],...
    'fontsize',12,...
    'string','Next >',...
    'callback',{@go_next, Seg});



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




set(S.ls, 'callback',{@list_call});

S.testFun = @testFun;
end

function [] = radio_call(varargin)
global S;
%S = varargin{3};  % Get the structure.

% Instead of switch, we could use num2str on:
% find(get(S.bg,'selectedobject')==S.rd)      (or similar)
% Note the use of findobj.  This is because of a BUG in MATLAB, whereby if
% the user selects the same button twice, the selectedobject property will
% not work correctly.
switch get(findobj(get(S.bg,'selectedobject')), 'String')
    case get(S.rd(1), 'String')
        disp('segments');
        flip_view(1);
    case get(S.rd(2), 'String')
        disp('rgb');
        flip_view(2);
    otherwise
        set(S.ed,'string','None!') % Very unlikely I think.
end
end

function [] = go_next(varargin)
global dataIdx;
global dataSize;

if dataIdx < dataSize
    dataIdx = dataIdx+1;
    updateFigure();
end
end

function [] = go_prev(varargin)
global dataIdx;


if dataIdx > 1
    dataIdx = dataIdx-1;
    updateFigure();
end
end

function [] = initFigure()
    global dataIdx;
    global plyData;
    global segData;
    global S;
    global dataRange;
    
    for i=1:4
        V2 = plyData.v{i};%coloredCloud(1:3, :);
        C2 = plyData.c{i};%coloredCloud(4:6, :);
        SV2 = segData{i}.SegVertex;
        SC2 = segData{i}.SegColor;
        Seg2 = segData{i};
        S.axSubt{i} = {};
        for j=1:length(Seg2.group)
            if j==1
                hold(S.axSub(i), 'off');
            end

            S.axSubt{i}{j} = scatter3(S.axSub(i), V2(1,Seg2.group{j})', V2(2,Seg2.group{j})', ...
                V2(3,Seg2.group{j})', 8, C2(:,Seg2.group{j})', 'filled');
            if j==1
                hold(S.axSub(i), 'on');
            end
            set(S.axSubt{i}{j}, 'Visible', 'off');


            S.axSubt_seg{i}{j} = scatter3(S.axSub(i), SV2(1,Seg2.group{j})', SV2(2,Seg2.group{j})', ...
                SV2(3,Seg2.group{j})', 8, SC2(:, Seg2.group{j})', 'filled');
            set(S.axSubt_seg{i}{j}, 'Visible', 'off');
        end

        view(S.axSub(i), [245 65]);
        grid(S.axSub(i), 'on');
        axis(S.axSub(i), 'equal');

        xlim(S.axSub(i), [dataRange.min(1) dataRange.max(1)]);
        ylim(S.axSub(i), [dataRange.min(2) dataRange.max(2)]);
        zlim(S.axSub(i), [dataRange.min(3) dataRange.max(3)]);
    end
end

function [] = updateFigure()
global dataIdx;
global plyData;
global segData;
global S;
global dataRange;

V = plyData.v{dataIdx};%coloredCloud(1:3, :);
C = plyData.c{dataIdx};%coloredCloud(4:6, :);
SV = segData{dataIdx}.SegVertex;
SC = segData{dataIdx}.SegColor;
Seg = segData{dataIdx};
%     global seg;
%     Seg = seg;
%     global S;

for i=1:length(Seg.group)
    if i==1
        hold(S.ax, 'off');
    end
    S.axt{i} = scatter3(S.ax, V(1,Seg.group{i})', V(2,Seg.group{i})', ...
        V(3,Seg.group{i})', 8, C(:,Seg.group{i})', 'filled');
    if i==1
        hold(S.ax, 'on');
    end
    S.axt_seg{i} = scatter3(S.ax, SV(1,Seg.group{i})', SV(2,Seg.group{i})', ...
        SV(3,Seg.group{i})', 8, SC(:, Seg.group{i})', 'filled');
    set(S.axt_seg{i}, 'Visible', 'off');
    
    if i==1
        hold(S.ax2, 'off');
    end
    S.axt2{i} = scatter3(S.ax2, V(1,Seg.group{i})', V(2,Seg.group{i})', ...
        V(3,Seg.group{i})', 8, C(:,Seg.group{i})', 'filled');
    if i==1
        hold(S.ax2, 'on');
    end
    set(S.axt2{i}, 'Visible', 'off');
    S.axt_seg2{i} = scatter3(S.ax2, SV(1,Seg.group{i})', SV(2,Seg.group{i})', ...
        SV(3,Seg.group{i})', 8, SC(:, Seg.group{i})', 'filled');
    set(S.axt_seg2{i}, 'Visible', 'off');
    
    %         visibleFigIdx2{i} = zeros(length(Vertex), 1);
end



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

view(S.ax, [245 65]);
view(S.ax2, [245 65]);
grid(S.ax, 'on');
grid(S.ax2, 'on');
axis(S.ax, 'equal');
axis(S.ax2, 'equal');
xlim(S.ax, [dataRange.min(1) dataRange.max(1)]);
ylim(S.ax, [dataRange.min(2) dataRange.max(2)]);
zlim(S.ax, [dataRange.min(3) dataRange.max(3)]);
xlim(S.ax2, [dataRange.min(1) dataRange.max(1)]);
ylim(S.ax2, [dataRange.min(2) dataRange.max(2)]);
zlim(S.ax2, [dataRange.min(3) dataRange.max(3)]);
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


function [] = flip_view(type)

global annotationResult;
global visibleGroupLeftIdx;
global selectedObjectIndex;
global S;
global dataIdx;

if type == 1
    for i=1:length(visibleGroupLeftIdx{dataIdx})
        if visibleGroupLeftIdx{dataIdx}(i) == 1
            set(S.axt{i}, 'Visible', 'on');
            set(S.axt_seg{i}, 'Visible', 'off');
        end
    end
    
    visibleSegGroup = annotationResult{dataIdx}{selectedObjectIndex}.group;
    for i = 1:length(visibleSegGroup)
        if visibleSegGroup(i) == 1
            set(S.axt2{i}, 'Visible', 'on');
            set(S.axt_seg2{i}, 'Visible', 'off');
        end
    end
elseif type == 2
    for i=1:length(visibleGroupLeftIdx{dataIdx})
        if visibleGroupLeftIdx{dataIdx}(i) == 1
            set(S.axt{i}, 'Visible', 'off');
            set(S.axt_seg{i}, 'Visible', 'on');
        end
    end
    
    visibleSegGroup = annotationResult{dataIdx}{selectedObjectIndex}.group;
    for i = 1:length(visibleSegGroup)
        if visibleSegGroup(i) == 1
            set(S.axt2{i}, 'Visible', 'off');
            set(S.axt_seg2{i}, 'Visible', 'on');
        end
    end
    
elseif type == 3
    
end

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
global dataIdx;

disp(sprintf('updateLeftFig Called!, prev:%d, cur:%d', prevSelectedObjectIndex, selectedObjectIndex));

for i=1:length(visibleGroupLeftIdx{dataIdx})
    if visibleGroupLeftIdx{dataIdx}(i) == 1
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
global dataIdx;

disp(sprintf('updateRightFig Called!, prev:%d, cur:%d', prevSelectedObjectIndex, selectedObjectIndex));

if prevSelectedObjectIndex ~= -1
    prevVisibleSegGroup = annotationResult{dataIdx}{prevSelectedObjectIndex}.group;
    for i = 1:length(prevVisibleSegGroup)
        if prevVisibleSegGroup(i) == 1
            set(S.axt2{i}, 'Visible', 'off');
            set(S.axt_seg2{i}, 'Visible', 'off');
            
        end
    end
end

visibleSegGroup = annotationResult{dataIdx}{selectedObjectIndex}.group;
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


