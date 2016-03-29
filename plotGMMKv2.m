function [h] = plotGMMKv2( params )
%PLOTGMMK Summary of this function goes here
%   Detailed explanation goes here
%     setId=1;
    %GMM group index
%     gmmk = 1;

    setId = params.setId;
    gmmk = params.gmmk;
    assignments = params.assignments;
    maxIdx = params.maxIdx;
    fv = params.fv;
    I2 = params.I2;
    
    

    if ~isfield(params, 'type')
        type = 1;
    else
        type = params.type;
    end
    
    if isfield(params, 'gmmk_gmmidx') && ~isempty(params.gmmk_gmmidx)
        gmmk_gmmidx = params.gmmk_gmmidx;
    else
        
        gmmk_gmmidx = [];
        for gmk = gmmk
        %included GMM index in goup GMMK
        gmmk_gmmidx = [gmmk_gmmidx find(assignments==gmk)];
        end
    end
    
       
   
    %included sample index
    gmmk_sampleIdx = find(ismember(maxIdx{setId}, gmmk_gmmidx));
    %included original Point cloud index
    gmmk_Idx = find(ismember(I2{setId}, gmmk_sampleIdx));

    fv2.Faces =  fv{setId}.Faces;
    fv2.Vertices = fv{setId}.Vertices;
    
    % 
    h = figure;
    
    %Color type
    if type == 1
        fv2.FaceVertexCData = fv{setId}.FaceVertexCData;
        if isempty(fv2.Faces)
            interval = 10;
            fv2.FaceVertexCData(gmmk_Idx,:) = bsxfun(@plus,fv2.FaceVertexCData(gmmk_Idx,:),[0.5 0 0]);
            scatter3(fv2.Vertices(1:interval:end,1), fv2.Vertices(1:interval:end,2), ...
                fv2.Vertices(1:interval:end,3), 8, fv2.FaceVertexCData(1:interval:end, :), 'filled');
            view(28, -20);
            axis equal;
        else
            fv2.FaceVertexCData = fv{setId}.FaceVertexCData;%sColor(:,I2{id})';
            fv2.FaceVertexCData(gmmk_Idx,:) = bsxfun(@plus,fv2.FaceVertexCData(gmmk_Idx,:),[0.5 0 0]);
            patch(fv2,'FaceColor','flat','EdgeColor','flat',...
                  'MarkerFaceColor','flat', ...
                  'FaceLighting',    'gouraud',     ...
                     'AmbientStrength', 0.15);
              title('color from Sampling');     
            axis vis3d;
        end
    %Dull type
    elseif type == 2
        if isempty(fv2.Faces)
            interval = 10;
            fv2.FaceVertexCData = zeros(size(fv{setId}.Vertices));%sColor(:,I2{id})';
            fv2.FaceVertexCData = bsxfun(@plus, fv2.FaceVertexCData, [0.7 0.7 0.9]);
            fv2.FaceVertexCData(gmmk_Idx,:) = bsxfun(@plus,bsxfun(@times, fv2.FaceVertexCData(gmmk_Idx,:), [0 0 0]),[0.9 0.4 0.4]);
            scatter3(fv2.Vertices(1:interval:end,1), fv2.Vertices(1:interval:end,2), ...
                fv2.Vertices(1:interval:end,3), 8, fv2.FaceVertexCData(1:interval:end, :), 'filled');
            view(28, -20);
            axis equal;
        else
            fv2.FaceVertexCData = zeros(size(fv{setId}.Vertices));%sColor(:,I2{id})';
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

end

