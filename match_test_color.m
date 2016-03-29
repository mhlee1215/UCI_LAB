% close;
% clear;

addpath(genpath('libs'));
run('libs2/vlfeat-0.9.20-bin/vlfeat-0.9.20/toolbox/vl_setup');


fv = {};
fprintf('data read...');
for data_idx = 1:2
    data_path = sprintf('data/desk_color2');
%     data_name = sprintf('MeshedReconstruction_%d.stl', data_idx);
%     crop_path = sprintf('%s/crop', data_path);
%     fv{data_idx} = stlread(sprintf('%s/%s', data_path, data_name));
    
    
    fprintf('loading.. %d/%d', data_idx, 1);
    [tri, pts, data, comments] = ply_read(sprintf('%s/MeshedReconstruction_%d.ply', data_path, data_idx), 'tri');
    
    fprintf('end\n');
    v = [data.vertex.x data.vertex.y data.vertex.z];
    f = tri';
    c = [data.vertex.red data.vertex.green data.vertex.blue];

    fv{data_idx}.faces = f;
    fv{data_idx}.vertices = v;
    fv{data_idx}.FaceVertexCData = c./255;
    
end
fprintf('fin \n');

fprintf('reducing..');
fv{1} = reducepatch(fv{1}, .2);
fv{2} = reducepatch(fv{2}, .2);
fprintf('end\n');

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




% DescrData = pdist2(img1_desc_types,img2_desc_types);
% 
% [aa,NN_Data]=sort(DescrData,2,'ascend');
% NN_Data=NN_Data(:,1:NN_Number);

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

% x = [vertex11(pts1,1) vertex22(pts2(NN_Data(:,1)), 1)]';
% y = [vertex11(pts1,2) vertex22(pts2(NN_Data(:,1)), 2)]';
% z = [vertex11(pts1,3) vertex22(pts2(NN_Data(:,1)), 3)]';
% 
% % validFlag = (aa(:,1) > 0);
% 
% validIdx = find((aa(:,1) > 0).* (aa(:,1) ./ aa(:,2) < 0.7) .* (aa(:,1) < 0.1));
% plot3(x(:,validIdx), y(:,validIdx), z(:,validIdx));



% return;

% [matches, scores] = vl_ubcmatch(img1_desc_types(1:100,:)',img2_desc_types(1:100,:)', 2) ;

[matches, scores] = vl_ubcmatch(img1_desc_types',img2_desc_types', 2) ;
numMatches = size(matches,2) ;

% x = [vertex11(pts1(matches(1,:)),1) vertex22(pts2(matches(2,:)),1)]';
% y = [vertex11(pts1(matches(1,:)),2) vertex22(pts2(matches(2,:)),2)]';
% z = [vertex11(pts1(matches(1,:)),3) vertex22(pts2(matches(2,:)),3)]';
% plot3(x(:,:), y(:,:), z(:,:));
% 
% return;


X1 = vertex11(matches(1,:), :)';%f1(1:2,matches(1,:)) ; X1(3,:) = 1 ;
X2 = vertex11(matches(2,:), :)';%f2(1:2,matches(2,:)) ; X2(3,:) = 1 ;

% --------------------------------------------------------------------
%                                         RANSAC with homography model
% --------------------------------------------------------------------
clear H score ok ;
for t = 1:500
  % estimate homograpyh
  subset = vl_colsubset(1:numMatches, 4) ;
  A = [] ;
  for i = subset
    A = cat(1, A, kron(X1(:,i)', vl_hat(X2(:,i)))) ;
  end
  [U,S,V] = svd(A) ;
  H{t} = reshape(V(:,9),3,3) ;

  % score homography
  X2_ = H{t} * X1 ;
  du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
  dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
  ok{t} = sqrt(du.*du + dv.*dv) < 0.05 ;
  scoreA(t) = sum(ok{t}) ;
end

[score, best] = max(scoreA) ;
H = H{best} ;
ok = ok{best} ;



% x = [vertex11(pts1(matches(1,ok)),1) vertex22(pts2(matches(2,ok)),1)]';
% y = [vertex11(pts1(matches(1,ok)),2) vertex22(pts2(matches(2,ok)),2)]';
% z = [vertex11(pts1(matches(1,ok)),3) vertex22(pts2(matches(2,ok)),3)]';
x = [vertex11(pts1(matches(1,:)),1) vertex22(pts2(matches(2,:)),1)]';
y = [vertex11(pts1(matches(1,:)),2) vertex22(pts2(matches(2,:)),2)]';
z = [vertex11(pts1(matches(1,:)),3) vertex22(pts2(matches(2,:)),3)]';
plot3(x(:,:), y(:,:), z(:,:));

