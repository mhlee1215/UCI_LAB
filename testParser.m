function [ p ] = testParser( input1, varargin )
%TESTPARSER Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
% p.addOptional('ntimes',1);

p.addParameter('binCapacity2',-1);
p.addParameter('binCapacity',-1);

p.parse(varargin{:})



end

