function [ mats ] = read4x4Mats( path )
%READ4X4MATS Summary of this function goes here
%   Detailed explanation goes here
    fileID = fopen(path, 'r');
    mat = fscanf(fileID, '%f');
    mat = reshape(mat, 4, length(mat)/4)';    
    
    mats = {};
    for pi = 1:size(mat, 1)/4
        M.R = mat(4*(pi-1)+1:4*(pi-1)+3, 1:3);
        M.t = mat(4*(pi-1)+1:4*(pi-1)+3, 4);
        mats{end+1} = M;
    end
end

