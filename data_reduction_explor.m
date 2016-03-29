
addpath(genpath('libs'));
d_min = 0.4;

center_roi = {};
center_roi{1} = [-0.02734 0.07161 -0.5351];
center_roi{2} = [-0.04818 0.01406 -0.5128];
center_roi{3} = [0.02474 0.07585 -0.4815];
center_roi{4} = [-0.03516 0.08076 -0.4893];

for data_idx = 5


data_path = sprintf('data/desk/%d', data_idx);
data_name = 'MeshedReconstruction.stl';
crop_path = sprintf('%s/crop', data_path);

result_path = sprintf('%s/%s', crop_path, data_name);
fv = stlread(result_path);

figure;
patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
view([-135 35]);

% figure;
% fv2 = reducepatch(fv, 0.9);
% nfv = reducepatch(fv.faces, fv.vertices, 0.01);
% patch(fv2,'FaceColor',       [0.8 0.8 1.0], ...
%          'EdgeColor',       'none',        ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);
% 
% % Add a camera light, and tone down the specular highlighting
% camlight('headlight');
% material('dull');
% 
% % Fix the axes scaling, and set a nice view angle
% axis('image');
% view([-135 35]);


end


