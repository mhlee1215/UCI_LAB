% close;
% clear;

addpath(genpath('libs'));
run('libs2/vlfeat-0.9.20-bin/vlfeat-0.9.20/toolbox/vl_setup');


fv = {};
fprintf('data read...');
for data_idx = 1:2
    data_path = sprintf('data/desk/%d', data_idx);
    data_name = 'MeshedReconstruction.stl';
    crop_path = sprintf('%s/crop', data_path);
    fv{data_idx} = stlread(sprintf('%s/%s', crop_path, data_name));
end
fprintf('fin \n');

fv{1} = reducepatch(fv{1}, 1);
fv{2} = reducepatch(fv{2}, 1);

vertex1 = fv{1}.vertices;
faces1 = fv{1}.faces;
vertex2 = fv{2}.vertices;
faces2 = fv{2}.faces;


bin = 3;
scale_rat = 1;
NN_Number = 5;

% [ img1_desc_types, detectedPts1, locations1,pts1 ] = create_mesh_sift_features(vertex1, faces1, scale_rat);
tic;
[ img1_desc_types, ~,~,pts1 ] = create_mesh_spin_features(vertex1, faces1,bin,scale_rat);
toc;

% [ img2_desc_types, detectedPts2, locations2,pts2] = create_mesh_sift_features(vertex2, faces2, scale_rat);
[ img2_desc_types, ~, ~,pts2] = create_mesh_spin_features(vertex2, faces2,bin,scale_rat);

DescrData = pdist2(img1_desc_types,img2_desc_types);

[aa,NN_Data]=sort(DescrData,2,'ascend');
NN_Data=NN_Data(:,1:NN_Number);

% return;

dist_between = mean(max(vertex1) - min(vertex1) )*2;
vertex11 = vertex1;
vertex11(:,3) = vertex11(:,3).*1;
vertex22 = vertex2+repmat([dist_between 0 0], size(vertex2, 1), 1);
vertex22(:,3) = vertex22(:,3).*1;

fv1.vertices = [vertex11 ; vertex22];
fv1.faces = [faces1 ; faces2+size(vertex11, 1)];
figure;
patch(fv1,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
view([-135 35]);

hold on;
scatter3(vertex11(pts1,1), vertex11(pts1,2), vertex11(pts1,3), '+');
scatter3(vertex22(pts2,1), vertex22(pts2,2), vertex22(pts2,3), '*');

% x = [0 , 3; -1, -5]';
% y = [0 , 3; -1, -5]';
% z = [0 , 3; -1, -5]';

x = [vertex11(pts1,1) vertex22(pts2(NN_Data(:,1)), 1)]';
y = [vertex11(pts1,2) vertex22(pts2(NN_Data(:,1)), 2)]';
z = [vertex11(pts1,3) vertex22(pts2(NN_Data(:,1)), 3)]';

validIdx = find((aa(:,1) > 0).* (aa(:,1) ./ aa(:,2) < 0.7) .* (aa(:,1) < 0.03));

plot3(x(:,validIdx), y(:,validIdx), z(:,validIdx));


