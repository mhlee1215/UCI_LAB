function myCallbackClickA3DPointRight(src, eventData, pointCloud, C, Seg)
% CALLBACKCLICK3DPOINT mouse click callback function for CLICKA3DPOINT
%
%   The transformation between the viewing frame and the point cloud frame
%   is calculated using the camera viewing direction and the 'up' vector.
%   Then, the point cloud is transformed into the viewing frame. Finally,
%   the z coordinate in this frame is ignored and the x and y coordinates
%   of all the points are compared with the mouse click location and the 
%   closest point is selected.
%
%   Babak Taati - May 4, 2005
%   revised Oct 31, 2007
%   revised Jun 3, 2008
%   revised May 19, 2009

disp('here is Right');

global visibleGroupLeftIdx;
global selectedObjectIndex;
global annotationResult;
global prevGroup;
global visibleFigIdx;
global visibleFigIdx2;
global S;
global dataIdx;

visibleFigIdxRight = visibleFigIdx2{dataIdx}{selectedObjectIndex};
hh = S.ax2;

point = get(hh, 'CurrentPoint'); % mouse click position

camPos = get(hh, 'CameraPosition'); % camera position
camTgt = get(hh, 'CameraTarget'); % where the camera is pointing to

camDir = camPos - camTgt; % camera direction
camUpVect = get(hh, 'CameraUpVector'); % camera 'up' vector

% build an orthonormal frame based on the viewing direction and the 
% up vector (the "view frame")
zAxis = camDir/norm(camDir);    
upAxis = camUpVect/norm(camUpVect); 
xAxis = cross(upAxis, zAxis);
yAxis = cross(zAxis, xAxis);

rot = [xAxis; yAxis; zAxis]; % view rotation 

% the point cloud represented in the view frame
disp(sprintf('Visible num : %d', sum(visibleFigIdxRight)));

visibleIndex = find(visibleFigIdxRight);
rotatedPointCloud = rot * pointCloud(:,visibleIndex); 

% the clicked point represented in the view frame
rotatedPointFront = rot * point' ;

% find the nearest neighbour to the clicked point 
[pointCloudIndex, distance] = dsearchn(rotatedPointCloud(1:2,:)', ... 
    rotatedPointFront(1:2));

h = findobj(hh,'Tag','pt'); % try to find the old point

if(max(max(abs(point))) > 30 || distance > 0.1) 
%     disp(sprintf('hihi, %f, %f', max(max(abs(point))), distance));
    if ~isempty(h)
        delete(h); % delete the previously selected point
    end
    prevGroup = 0;
    return; 
end


% prevGroup = findobj(hh,'Tag','prevGroup'); % try to find the old point
% prevGroup
selectedPoint = pointCloud(:, visibleIndex(pointCloudIndex)); 

selectedGroupNumber = Seg.idx(visibleIndex(pointCloudIndex));
disp(sprintf('group number : %d selected', selectedGroupNumber));
selectedGroupIdx = Seg.group{selectedGroupNumber};
selectedGroupPoint = pointCloud(:, selectedGroupIdx);

if isempty(h) % if it's the first click (i.e. no previous point to delete)
   disp('111hihi'); 
    % highlight the selected point
%     h = plot3(selectedPoint(1,:), selectedPoint(2,:), ...
%         selectedPoint(3,:), 'r.', 'MarkerSize', 20); 
    
    h = plot3(selectedGroupPoint(1,:), selectedGroupPoint(2,:), ...
        selectedGroupPoint(3,:), 'r.', 'MarkerSize', 20); 
    
    set(h,'Tag','pt'); % set its Tag property for later use   
    set(h, 'ButtonDownFcn', {@myCallbackClickA3DPointRight, pointCloud, C, Seg}); 
    prevGroup = selectedGroupNumber;
else % if it is not the first click
disp('2hihi');
    delete(h); % delete the previously selected point
%     prevGroup = 0;
    % highlight the newly selected point
%     h = plot3(selectedPoint(1,:), selectedPoint(2,:), ...
%         selectedPoint(3,:), 'r.', 'MarkerSize', 20);  

     if prevGroup ~= selectedGroupNumber 
         h = plot3(selectedGroupPoint(1,:), selectedGroupPoint(2,:), ...
            selectedGroupPoint(3,:), 'r.', 'MarkerSize', 20); 
         set(h,'Tag','pt');  % set its Tag property for later use
         set(h, 'ButtonDownFcn', {@myCallbackClickA3DPointRight, pointCloud, C, Seg}); 
         prevGroup = selectedGroupNumber;
     else
         prevGroup = 0;
         
%          set(S.axt{selectedGroupNumber}, 'Visible', 'off');
%          set(S.axt_seg{selectedGroupNumber}, 'Visible', 'off');
         
         set(S.axt2{selectedGroupNumber}, 'Visible', 'off');
         set(S.axt_seg2{selectedGroupNumber}, 'Visible', 'off');
         visibleFigIdx{dataIdx}(selectedGroupIdx) = 1;
         visibleGroupLeftIdx{dataIdx}(selectedGroupNumber) = 1;
         
         set(S.axSubt{dataIdx}{selectedGroupNumber}, 'Visible', 'off');
         set(S.axSubt_seg{dataIdx}{selectedGroupNumber}, 'Visible', 'off');
         
         disp(sprintf('group number : %d deleted', selectedGroupNumber));
         annotationResult{dataIdx}{selectedObjectIndex}.group(selectedGroupNumber) = 0;
         %Update fig... TODO!
%          S.axt2{selectedGroupNumber} = scatter3(S.ax2, pointCloud(1,selectedGroupIdx), pointCloud(2,selectedGroupIdx), ...
%             pointCloud(3,selectedGroupIdx), 8, C(:,selectedGroupIdx)', 'filled'); hold on;
%         
         S.updateLeftFig();
     end
end



pointCloud = [];
% set(selectedGroupNumber, 'Tag', 'prevGroup');
fprintf('you clicked on point number %d\n', pointCloudIndex);
