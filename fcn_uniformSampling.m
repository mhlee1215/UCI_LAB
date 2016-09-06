function [ vs, cs, ns, sIndex ] = fcn_uniformSampling( v, c, n, density )
%FCN_UNIFORMSAMPLING Summary of this function goes here
%   Detailed explanation goes here

        [vs, sIndex, cs] = uniformSubSample(v', density, c');
        [~, ~, ns] = uniformSubSample(v', density, n');

        vs = vs';
        cs = cs';
        ns = ns';
        
end

