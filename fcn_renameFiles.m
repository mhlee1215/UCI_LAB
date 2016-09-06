function [ ] = fcn_renameFiles( root, identifier, src, dst )
%FCN_RENAMEFILES Summary of this function goes here
%   Detailed explanation goes here
% fcn_renameFiles('/home/mhlee/data_from_odroid/complete', 'LAB_1-2016-07-29_12_39.klg', 'LAB_1', 'LAB_2');

% Get all PDF files in the current folder
files = dir(sprintf('%s/*%s*', root, identifier));
% Loop through each
for id = 1:length(files)
        oldStr = files(id).name;
        newStr = strrep(files(id).name, src, dst);
        disp(sprintf('string %s -> %s', oldStr, newStr));
        oldPath = sprintf('%s/%s', root, oldStr);
        newPath = sprintf('%s/%s', root, newStr);
        movefile(oldPath, newPath);
end


end

