function [ vertices_sub ] = randSubSample( vertices, sampleNum )
%SUBSAMPLE Summary of this function goes here
%   Detailed explanation goes here
rand_idx = randperm(size(vertices, 1));
vertices_sub = vertices(rand_idx(1:sampleNum), :)';

end

