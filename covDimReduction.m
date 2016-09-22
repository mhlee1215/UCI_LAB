function [ Q2 ] = covDimReduction( Q )
%COVDIMREDUCTION Summary of this function goes here
%   Detailed explanation goes here

[U, S, V] = svd(Q);
Q2 = U*[S(1,1) 0 0 ; 0 S(2,2) 0 ; 0 0 S(3, 3)*0]*V';


end

