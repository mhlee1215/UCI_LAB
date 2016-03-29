function [ A, maxIdx, assignments, GMM_mean_color, A_binary, A_binary_c, isStatic ] = posteriorAnalysis( a, X, C2, K, clusterSize, var_threshold )
%POSTERIORANALYSIS Summary of this function goes here
%   Detailed explanation goes here

d_size = length(a);

A = zeros(K, length(a));

%A2 = zeros(K, d_size);
%A3 = zeros(K, d_size);
%A4 = zeros(K, d_size);%zeros(K, 1);

maxIdx = {};

for ii=1:length(a)
    [m, i] = max(a{ii}');
    maxIdx{ii} = i;
    gmmCount = accumarray(i', 1);
    A(1:length(gmmCount), ii) = gmmCount;

%     [sortedValues,sortIndex] = sort(a{ii}','descend'); 


%     firstIdx = sortIndex(1, :)';
%     firstVal = sortedValues(1, :)';
%     secondIdx = sortIndex(2, :)';
%     secondVal = sortedValues(2, :)';

    % KKK = accumarray(firstIdx, firstVal);
    % KKK2 = accumarray(secondIdx, secondVal);
    % 
    % A4 = [A4 KKK+KKK2];
end
% AA = bsxfun(@rdivide, A, max(A')');
% figure; plot(AA');
% h=figure; plot(A'); title(sprintf('Group size of each GMM (#K : %d)', K));
% saveas(h, sprintf('%s/groupassign_k%d', resultPath, K));
% close(h);
% figure; plot(A2');
% figure; plot(A3');
% figure; plot(A4');

% gmm_score = AA;
% 
% cluster = [];
% for i=1:K
%     cluster = [cluster vl_kmeans(gmm_score(i,:), 2)];
% end
% figure; plot(cluster)


% return;

%rGMM_color = bsxfun(@plus, rand(K, 3)./2, [0.3 0.3 0.3]);
% for id=1:d_size
%     fv3.Faces = fv{id}.Faces;
%     fv3.Vertices = fv{id}.Vertices;
% 
%     sColor = C2{id};
% 
%     mIdx = maxIdx{id};
%     sRMean2 = accumarray(mIdx', sColor(1,:)', [], @mean);
%     sGMean2 = accumarray(mIdx', sColor(2,:)', [], @mean);
%     sBMean2 = accumarray(mIdx', sColor(3,:)', [], @mean);
% 
%     GMM_color = [sRMean2 sGMean2 sBMean2];
%     sampleColor = GMM_color(mIdx, :)';
% 
%     sIndex = I2{id};
%     fv3.FaceVertexCData = sampleColor(:,sIndex)';
% 
%     h = figure;
%     patch(fv3,'FaceColor','flat','EdgeColor','flat',...
%           'MarkerFaceColor','flat', ...
%           'FaceLighting',    'gouraud',     ...
%              'AmbientStrength', 0.15);
%     title('color from GMM');     
%     axis vis3d;   
%     saveas(h, sprintf('%s/color_GMM_%d', resultPath, id));
%     close(h);
% 
%     
%     
%     
%     sampleColor = rGMM_color(mIdx, :)';
%     sIndex = I2{id};
%     fv3.FaceVertexCData = sampleColor(:,sIndex)';
% 
%     h = figure;
%     patch(fv3,'FaceColor','flat','EdgeColor','flat',...
%              'EdgeColor',       'none',        ...
%              'FaceLighting',    'gouraud',     ...
%              'AmbientStrength', 0.15);
%     title('color from GMM');     
%     camlight('headlight');
%     material('dull');
%     axis vis3d;   
%     saveas(h, sprintf('%s/rcolor_GMM_%d', resultPath, id));
%     close(h);
%     
% 
%     fv2.Faces =  fv{id}.Faces;
%     fv2.Vertices = fv{id}.Vertices;
%     sColor = C2{id};
%     fv2.FaceVertexCData = sColor(:,I2{id})';
%     % 
%     h = figure;
%     patch(fv2,'FaceColor','flat','EdgeColor','flat',...
%           'MarkerFaceColor','flat', ...
%           'FaceLighting',    'gouraud',     ...
%              'AmbientStrength', 0.15);
%       title('color from Sampling');     
%     axis vis3d;
%     saveas(h, sprintf('%s/color_sampling_%d.fig', resultPath, id));
%     close(h);
% end

GMM_sum_color = zeros(K, 3);
% GMM_mean_color = zeros(K, 3);
GMM_sum_index = zeros(K, 1);
for id=1:d_size
    sColor = C2{id};

    mIdx = maxIdx{id};
    sRMean2 = accumarray(mIdx', sColor(1,:)', [], @mean);
    sGMean2 = accumarray(mIdx', sColor(2,:)', [], @mean);
    sBMean2 = accumarray(mIdx', sColor(3,:)', [], @mean);

    GMM_color = [sRMean2 sGMean2 sBMean2];
    GMM_color_container = zeros(K, 3);
    GMM_color_container(1:size(GMM_color, 1), :) = GMM_color;
    
    GMM_sum_color = GMM_sum_color + GMM_color_container;
    GMM_sum_index = GMM_sum_index + (sum(GMM_color_container==0, 2) < 3);
end
GMM_mean_color = bsxfun(@rdivide, GMM_sum_color, GMM_sum_index);
GMM_nan_idx = isnan(GMM_mean_color);
GMM_mean_color(find(GMM_nan_idx)) = 0;


AX = [A' ;X.*sqrt(size(A,2)/3).*mean(mean(A)).*1.2 ; GMM_mean_color'.*sqrt(size(A,2)/3).*mean(mean(A))];
labels = '_AX3';
% [centers, assignments] = vl_kmeans(A', 16);

[~, assignments] = vl_kmeans(AX, clusterSize);

A2 = A';
A_binary = zeros(size(A2));

rmax = 5;
cmax = 5;

% h=figure('units','normalized','outerposition',[0 0 1 1]);
% pIdx = 1;
% ha = tight_subplot(rmax, cmax,[.03 .01],[.01 .05],[.01 .01]); 
% for rr=1:rmax
%     for cc=1:cmax
%         subplot(rmax, cmax, pIdx);   
%         
% %         axes(ha(pIdx));
%         hold on;
%         assinIdx = find(assignments == pIdx);
%         A2sub = A2(:, assinIdx);
%         [m, i] = max(mean(A2sub'));
%         plot(A2(:, find(assignments == pIdx)));
%         plot(mean(A2sub'), 'k', 'lineWidth', 2);
%         ylim([0 60]);
%         title(sprintf('%d-th, size : %d, var:%.2f', pIdx, sum(assignments == pIdx), var(mean(A2sub'))));
%         hold off;
%         pIdx = pIdx+1;
%     end
% end
% saveas(h, sprintf('%s/groupassign_kmeans%s', resultPath, labels));
% close(h);

% h=figure('units','normalized','outerposition',[0 0 1 1]);
% pIdx = 1;
% % ha = tight_subplot(rmax, cmax,[.03 .01],[.01 .05],[.01 .01]); 
% for rr=1:rmax
%     for cc=1:cmax
%         subplot(rmax, cmax, pIdx);   
%         
% %         axes(ha(pIdx));
%         hold on;
%         assinIdx = find(assignments == pIdx);
%         A2sub = A2(:, assinIdx);
%         [m, i] = max(mean(A2sub'));
%         plot(A2(:, find(assignments == pIdx)));
%         plot(mean(A2sub'), 'k', 'lineWidth', 2);
%         plot(mean(mean(A2sub')).*ones(size(mean(A2sub'))), 'g', 'lineWidth', 2);
%         ylim([0 60]);
%         title(sprintf('%d-th, size : %d, var:%.2f', pIdx, sum(assignments == pIdx), var(mean(A2sub'))));
%         hold off;
%         pIdx = pIdx+1;
%     end
% end
% saveas(h, sprintf('%s/groupassign_kmeans_threshold%s', resultPath, labels));
% close(h);

% h=figure('units','normalized','outerposition',[0 0 1 1]);

% var_threshold = 15;
% ha = tight_subplot(rmax, cmax,[.03 .01],[.01 .05],[.01 .01]); 

isStatic = [];
for pIdx = 1:clusterSize
    assinIdx = find(assignments == pIdx);
    A2sub = A2(:, assinIdx);
    if var(mean(A2sub')) > var_threshold
        %Binarize
        A_binary(:, assinIdx) = repmat((mean(A2sub') - mean(mean(A2sub')) > 0)', 1, length(assinIdx));
        isStatic(pIdx) = 0;
    else
        A_binary(:, assinIdx) = ones(size(A2sub));
        isStatic(pIdx) = 1;
    end
end

A_binary_c = {};
for i=1:size(A_binary, 1)
    A_binary_c{end+1} = A_binary(i, :);
end

end

