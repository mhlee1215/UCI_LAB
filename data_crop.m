
addpath(genpath('libs'));
d_min = 0.35;

center_roi = {};
center_roi{1} = [-0.02734 0.07161 -0.5351]; 
center_roi{2} = [-0.04818 0.01406 -0.5128];
center_roi{3} = [0.02474 0.07585 -0.4815]; 
center_roi{4} = [-0.03516 0.08076 -0.4893];

for data_idx = 1:4 
data_idx

data_path = sprintf('data/desk/%d', data_idx);
data_name = 'MeshedReconstruction';
data_ext = 'stl';
crop_path = sprintf('%s/crop', data_path);
mkdir(crop_path);
result_path = sprintf('%s/%s_%.2f.%s', crop_path, data_name, d_min, data_ext);

fv = stlread(sprintf('%s/%s.%s', data_path, data_name, data_ext));


% fv = stlread(result_path);
% figure;
% patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
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

% return;
c = center_roi{data_idx};
%c = [0.0664 0.7856 -0.5388];
% c = [0.09245 0.1343 -0.6091];
% c = [-0.02734 0.07161 -0.5351];
dist = fv.vertices - repmat(c, size(fv.vertices, 1), 1);
sq_dist = sqrt(sum(dist .^2, 2));
valid_idx = find(sq_dist < d_min);
v_crop = fv.vertices(valid_idx, :);

valid_face_idx = find(ismember(fv.faces(:,1), valid_idx).*ismember(fv.faces(:,2), valid_idx).*ismember(fv.faces(:,3), valid_idx));


k = valid_idx;
v = 1:length(valid_idx);
map = containers.Map(k, v);

f1_cvted = arrayfun(@(x) map(x), fv.faces(valid_face_idx, 1));
f2_cvted = arrayfun(@(x) map(x), fv.faces(valid_face_idx, 2));
f3_cvted = arrayfun(@(x) map(x), fv.faces(valid_face_idx, 3));

fv_crop.vertices = v_crop;
fv_crop.faces = [f1_cvted f2_cvted f3_cvted];

% figure;
% patch(fv_crop,'FaceColor',       [0.8 0.8 1.0], ...
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

stlwrite(result_path, fv_crop);
end


