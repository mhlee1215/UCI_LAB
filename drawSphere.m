function [ output_args ] = drawSphere( x, y, z, r )
%DRAWSPHERE Summary of this function goes here
%   Detailed explanation goes here

[xx,yy,zz] = sphere();
surf(xx*r+x, yy*r+y, zz*r+z);

end

