function [ v_crop, f_crop ] = cropFV( vertices, faces, valid_idx)
%CROPFV Summary of this function goes here
%   Detailed explanation goes here
    valid_face_idx = find(ismember(faces(:,1), valid_idx).*ismember(faces(:,2), valid_idx).*ismember(faces(:,3), valid_idx));
    v_crop = vertices(:,valid_idx);
    k = valid_idx;
    v = 1:length(valid_idx);
    map = containers.Map(k, v);

    f1_cvted = arrayfun(@(x) map(x), faces(valid_face_idx, 1));
    f2_cvted = arrayfun(@(x) map(x), faces(valid_face_idx, 2));
    f3_cvted = arrayfun(@(x) map(x), faces(valid_face_idx, 3));

    %fv_crop.vertices = v_crop;
    f_crop = [f1_cvted f2_cvted f3_cvted];

end

