clc;
close all;
clear all;

addpath(genpath('libs'));
load('match_exp_1.mat');


iidx = 1;
ii = IDX_SET(iidx);
ii_idx = find(IDX_1==ii);
ii_idx2 = find(IDX_2==ii);

v_m1 = [fv1_v(ii_idx, :)];
c_m1 = [fv1_c(ii_idx, :)];
% pclviewer([v_m1 c_m1]');

v_m2 = [fv2_v(ii_idx2, :)];
c_m2 = [fv2_c(ii_idx2, :)];
% pclviewer([v_m2 c_m2]');
v_m2_best = bsxfun(@plus, bestSubRSet{iidx}*v_m2', bestSubTSet{iidx})';
pclviewer([v_m1 c_m1 ; v_m2_best c_m2]');

% pclviewer([v_m1 c_m1 ; v_m2 c_m2]');
% 
% pclviewer([v_m1 c_m1]');
% pclviewer([v_m2 c_m2]');

v_m1_all = fv1_v;
c_m1_all = fv1_c;

v_m2_all = [];
c_m2_all = [];
for iidx=1%length(bestSubRSet)
    if isempty(bestSubRSet{iidx})
        continue;
    end
    
    ii = IDX_SET(iidx);
%     ii_idx = find(IDX_1==ii);
    ii_idx2 = find(IDX_2==ii);

%     v_m1 = [fv1_v(ii_idx, :)];
%     c_m1 = [fv1_c(ii_idx, :)];
    % pclviewer([v_m1 c_m1]');

    v_m2 = [fv2_v(ii_idx2, :)];
    c_m2 = [fv2_c(ii_idx2, :)];
    % pclviewer([v_m2 c_m2]');
    v_m2_best = bsxfun(@plus, bestSubRSet{iidx}*v_m2', bestSubTSet{iidx})';
    
    v_m2_all = [v_m2_all ; v_m2_best];
    c_m2_all = [c_m2_all ; c_m2];
    
end

pclviewer([v_m2_all c_m2_all]');