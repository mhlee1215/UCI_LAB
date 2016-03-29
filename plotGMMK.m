function [h] = plotGMMK( setId, gmmk, assignments, maxIdx, fv, I2, type )
%PLOTGMMK Summary of this function goes here
%   Detailed explanation goes here
%     setId=1;
    %GMM group index
%     gmmk = 1;


    if ~exist('type', 'var')
        type = 1;
    end
    
    gmmk_gmmidx = [];
    for gmk = gmmk
    %included GMM index in goup GMMK
    gmmk_gmmidx = [gmmk_gmmidx find(assignments==gmk)];
    end
    %included sample index
    gmmk_sampleIdx = find(ismember(maxIdx{setId}, gmmk_gmmidx));
    %included original Point cloud index
    gmmk_Idx = find(ismember(I2{setId}, gmmk_sampleIdx));

    fv2.Faces =  fv{setId}.Faces;
    fv2.Vertices = fv{setId}.Vertices;
    fv2.FaceVertexCData = fv{setId}.FaceVertexCData;%sColor(:,I2{id})';
    fv2.FaceVertexCData(gmmk_Idx,:) = bsxfun(@plus,fv2.FaceVertexCData(gmmk_Idx,:),[0.5 0 0]);
    % 
    h = figure;
    
    %Color type
    if type == 1
    
        patch(fv2,'FaceColor','flat','EdgeColor','flat',...
              'MarkerFaceColor','flat', ...
              'FaceLighting',    'gouraud',     ...
                 'AmbientStrength', 0.15);
          title('color from Sampling');     
        axis vis3d;
    
    %Dull type
    elseif type == 2
    
        fv2.FaceVertexCData = zeros(size(fv{setId}.FaceVertexCData));%sColor(:,I2{id})';
        fv2.FaceVertexCData = bsxfun(@plus, fv2.FaceVertexCData, [0.7 0.7 0.9]);
        fv2.FaceVertexCData(gmmk_Idx,:) = bsxfun(@plus,bsxfun(@times, fv2.FaceVertexCData(gmmk_Idx,:), [0 0 0]),[0.9 0.4 0.4]);
        patch(fv2, ...
             'FaceColor','flat','EdgeColor','flat', ...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);
        camlight('headlight');
        material('dull');
         axis vis3d;
    end

end

