function [ groups ] = Clustering( v, k, type )
%MYCLUSTERING Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3 || strcmp(type, 'hard')==1
    groups = ClusteringHard(v, k);
elseif strcmp(type, 'soft')==1
    groups = ClusteringSoft(v, k);
end

end

