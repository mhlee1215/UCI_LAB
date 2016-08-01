function [ mat ] = rotationMat3( theta, axis )
%ROTATIONMAT Summary of this function goes here
%   Detailed explanation goes here

mat = [];
if axis == 'x'
    mat = [1 0 0 ; 0 cos(theta) -sin(theta) ; 0 sin(theta) cos(theta)];
elseif axis == 'y'
    mat = [cos(theta) 0 sin(theta) ; 0 1 0 ; -sin(theta) 0 cos(theta)];
elseif axis == 'z'
    mat = [cos(theta) -sin(theta) 0 ; sin(theta) cos(theta) 0 ; 0 0 1];
end

end

