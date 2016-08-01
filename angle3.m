function [ theta ] = angle3( v1, v2 )
%ANGLE3 Summary of this function goes here
%   Detailed explanation goes here

theta = acos(dot(v1, v2) / (sqrt(sum(v1.^2))*sqrt(sum(v2.^2))));


end

