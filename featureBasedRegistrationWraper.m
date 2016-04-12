function [ bestSubR, bestSubT, d_seeds_12, f_12, d_seeds_22, f_22 ] = featureBasedRegistrationWraper( v_m1, c_m1, n_m1, v_m2, c_m2, n_m2, params )
%FEATUREBASEDREGISTRATIONWRAPER Summary of this function goes here
%   Detailed explanation goes here


    dataUnitSize = params.dataUnitSize;
    keypointUnitSize = params.keypointUnitSize;
    descriptorUnitSize = params.descriptorUnitSize;

    [v_m1_s2,~,c_m1_s2] = uniformSubSample(v_m1, dataUnitSize, c_m1);
    v_m1_s2 = v_m1_s2';
    c_m1_s2 = c_m1_s2';
    [v_m2_s2,~,c_m2_s2] = uniformSubSample(v_m2, dataUnitSize, c_m2);
    v_m2_s2 = v_m2_s2';
    c_m2_s2 = c_m2_s2';
    
    if ~isempty(n_m1)
        [~,~,n_m1_s2] = uniformSubSample(v_m1, dataUnitSize, n_m1);
        n_m1_s2 = n_m1_s2';
    else
        n_m1_s2 = [];
    end
    
    if ~isempty(n_m2)
        [~,~,n_m2_s2] = uniformSubSample(v_m2, dataUnitSize, n_m2);
        n_m2_s2 = n_m2_s2';
    else
        n_m2_s2 = [];
    end


    %Subsampling for dense descriptors
%     params_desc.normalRadius=0.2;
%     params_desc.searchRadius=0.8;
%     params_desc.searchK = 0;


    [v_m1_desc,~,c_m1_desc] = uniformSubSample(v_m1, descriptorUnitSize, c_m1);
    v_m1_desc = v_m1_desc';
    c_m1_desc = c_m1_desc';
    [v_m2_desc,~,c_m2_desc] = uniformSubSample(v_m2, descriptorUnitSize, c_m2);
    v_m2_desc = v_m2_desc';
    c_m2_desc = c_m2_desc';

    
    params_desc.normalRadius=params.normalRadius;
    params_desc.searchRadius=params.searchRadius;
    params_desc.searchK = params.searchK;

    if isempty(params.d_seeds_12)
        [d_seeds_12,~,~] = uniformSubSample(v_m1, keypointUnitSize, c_m1);
    else
        d_seeds_12 = params.d_seeds_12;
    end
    
    if isempty(params.f_12)
        f_12 = FPFHSDescriptor(v_m1_desc, c_m1_desc, d_seeds_12, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);
    else
        f_12 = params.f_12;
    end
    
%     figure; hold on; 
%     scatter3(v_m1_s2(:,1), v_m1_s2(:,2), v_m1_s2(:,3), 'k', 'filled'); 
%     scatter3(d_seeds_12(:,1), d_seeds_12(:,2), d_seeds_12(:,3), 'r', 'filled');
%     
%     scatter3(d_seeds_12(1,:), d_seeds_12(2,:), d_seeds_12(3,:), 'r', 'filled');
    
    % mean(sum(f_12==0, 1))
    if isempty(params.d_seeds_22)
        [d_seeds_22,~,~] = uniformSubSample(v_m2, keypointUnitSize, c_m2);
    else
        d_seeds_22 = params.d_seeds_22;
    end
    if isempty(params.f_22)
        f_22 = FPFHSDescriptor(v_m2_desc, c_m2_desc, d_seeds_22, params_desc.normalRadius, params_desc.searchRadius, params_desc.searchK);
    else
        f_22 = params.f_22;
    end
    % mean(sum(f_22==0, 1))

    if size(f_12, 2) < 3 || size(f_22, 2) < 3
        bestSubR = [];
        bestSubT = [];
        return;
    else
        [ bestSubR, bestSubT ] = featureBasedRegistration( v_m1_s2, c_m1_s2, n_m1_s2, v_m2_s2, c_m2_s2, n_m2_s2, d_seeds_12, f_12, d_seeds_22, f_22, 400 );    
    end
    
    


end

