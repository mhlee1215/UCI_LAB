function [ dataSetSampled ] = loadUniformSampling( dataSet, density )
%LOADUNIFORMSAMPLING Summary of this function goes here
%   Detailed explanation goes here

    dataSetSampled = {};
    progressbar2(0);
    progressbar2('Uniform sampling...');
    for i=1:length(dataSet)
        [vs, sIndex, cs] = uniformSubSample(dataSet{i}.v', density, dataSet{i}.c');
        [~, ~, ns] = uniformSubSample(dataSet{i}.v', density, dataSet{i}.n');
        
        dataSetSampled{i}.v = vs;
        dataSetSampled{i}.c = cs;
        dataSetSampled{i}.n = ns;
        dataSetSampled{i}.sIndex = sIndex;
        progressbar2(i/length(dataSet));
    end
end

