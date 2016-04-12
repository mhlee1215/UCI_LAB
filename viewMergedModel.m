function [ output_args ] = viewMergedModel( fv1_v, fv1_c, fv2_v, fv2_c, bestSubRSet, bestSubTSet, group_1, group_2 )
%VIEWMERGEDMODEL Summary of this function goes here
%   Detailed explanation goes here

v_m1_all = [];
c_m1_all = [];

v_m2_all = fv2_v;
c_m2_all = fv2_c;
for ii=1:length(group_1)
%     iidx
%     ii = IDX_SET(iidx);
    
    if ii > length(bestSubRSet) || isempty(bestSubRSet{ii})
        RR = eye(3);
        TT = [0 ; 0 ; 0];
%         continue;
    else
        RR = bestSubRSet{ii};
        TT = bestSubTSet{ii};
    end
    
    ii_idx = group_1{ii};%find(IDX_1==ii);
%     ii_idx2 = group_2{ii};%find(IDX_2==ii);

    v_m1 = [fv1_v(ii_idx, :)];
    c_m1 = [fv1_c(ii_idx, :)];
    % pclviewer([v_m1 c_m1]');

%     v_m2 = [fv2_v(ii_idx2, :)];
%     c_m2 = [fv2_c(ii_idx2, :)];
    % pclviewer([v_m2 c_m2]');
    v_m1_best = (RR'*bsxfun(@minus, v_m1', TT))';
    
    v_m1_all = [v_m1_all ; v_m1_best];
    c_m1_all = [c_m1_all ; c_m1];
    
end

pclviewer([v_m1_all c_m1_all]');

% pclviewer([v_m2_all c_m2_all]');
pclviewer([v_m1_all c_m1_all;v_m2_all c_m2_all]');

end

