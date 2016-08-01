function [ ] = drawMesh( v2, f2, c2 )
%DRAWMESH Summary of this function goes here
%   Detailed explanation goes here

    figure;
    fv1.vertices = v2;
    fv1.faces = f2;
    fv1.FaceVertexCData = c2;
    patch(fv1, ...
       'FaceColor','flat','EdgeColor','flat', ...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);
        camlight('headlight');
        material('dull');
    axis equal;

end

