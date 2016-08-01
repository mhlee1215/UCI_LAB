function [ m ] = RotMat( v1, v2 )
%ROTATIONMATMATLAB Summary of this function goes here
%   matlab in-built rotation matrix generation..

m = vrrotvec2mat(real(vrrotvec(v1, v2)));


end

