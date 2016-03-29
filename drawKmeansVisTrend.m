function [ h ] = drawKmeansVisTrend( A2, assignments, clusterSize )
%DRAWKMEANSVISTREND Summary of this function goes here
%   Detailed explanation goes here
rmax = (sqrt(clusterSize));
cmax = (sqrt(clusterSize));

h=figure('units','normalized','outerposition',[0 0 1 1]);
pIdx = 1;
% ha = tight_subplot(rmax, cmax,[.03 .01],[.01 .05],[.01 .01]); 
for rr=1:rmax
    for cc=1:cmax
        subplot(rmax, cmax, pIdx);   
        
%         axes(ha(pIdx));
        hold on;
        assinIdx = find(assignments == pIdx);
        A2sub = A2(:, assinIdx);
        [m, i] = max(mean(A2sub'));
        plot(A2(:, find(assignments == pIdx)));
        plot(mean(A2sub'), 'k', 'lineWidth', 2);
        plot(mean(mean(A2sub')).*ones(size(mean(A2sub'))), 'g', 'lineWidth', 2);
        ylim([0 60]);
        title(sprintf('%d-th, size : %d, var:%.2f', pIdx, sum(assignments == pIdx), var(mean(A2sub'))));
        hold off;
        pIdx = pIdx+1;
    end
end
% saveas(h, sprintf('%s/groupassign_kmeans%s', resultPath, labels));
% close(h);

end

