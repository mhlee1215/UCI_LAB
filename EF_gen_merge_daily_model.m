

clc;
close all;
% clear all;

addpath(genpath('libs'));

dataRoot = '/home/mhlee/data_from_odroid/complete/';
destPath = '/home/mhlee/data_from_odroid/merged/';


%Read All category set
fileListAll = dir(sprintf('%sLAB*cvt.ply', dataRoot));
category_set = {};
strCategory = '';
date_set = {};
data_set = {};
strDate = '';
dateMap = containers.Map();

for i=1:length(fileListAll)

    nameParts = strsplit(fileListAll(i).name, '-');
    curStrCategory = nameParts{1};
    
    nameParts2 = strsplit(nameParts{4}, '_');
    curStrDate = sprintf('%s-%s-%s', nameParts{2}, nameParts{3}, nameParts2{1});
        
    if ~isKey(dateMap, curStrDate)
        dateMap(curStrDate) = {};
    else
%         disp('hi');
    end
    
    curNameSet = dateMap(curStrDate);
    curNameSet{end+1} = fileListAll(i).name;
    dateMap(curStrDate) = curNameSet;    
end

colors = distinguishable_colors(30, [0.3 0.3 0.3]);
keySet = dateMap.keys;
for ki = 1:length(keySet)
    
    key = keySet{ki};
    fileList = dateMap(key);
    
    
    Rv = [];
    Rc = [];
    Rs = [];
    
    prevCatName = '';
    cnt = 1;
    for i=1:length(fileList)
        
        nameParts = strsplit(fileList{i}, '-');
        curStrCategory = nameParts{1};
               
        if strcmp(prevCatName, curStrCategory) == 1
            continue;
        end
        prevCatName = curStrCategory;
        
        model = ply_readWrap ( sprintf('%s%s', dataRoot, fileList{i}));
        Rv = [Rv ; model.v];
        Rc = [Rc ; model.c./255];
        Rs = [Rs ; repmat(colors(cnt, :), size(model.v, 1), 1)];
        
        cnt = cnt + 1;
    end
    
    
    % pclviewer(pcl_model);
    fileName = sprintf('%s_merged_rgb2', key);
    fcn_saveUniformSizeModel( Rv, Rc, [], destPath, fileName, 0 );
    fileName = sprintf('%s_merged_color2', key);
    fcn_saveUniformSizeModel( Rv, Rs, [], destPath, fileName, 0 );
end

