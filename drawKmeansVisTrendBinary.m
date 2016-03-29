function [ h ] = drawKmeansVisTrendBinary( A_binary, assignments, clusterSize )
%DRAWKMEANSVISTRENDBINARY Summary of this function goes here
%   Detailed explanation goes here

rmax = sqrt(clusterSize);
cmax = sqrt(clusterSize);

h=figure('units','normalized','outerposition',[0 0 1 1]);

pIdx = 1;
for rr=1:rmax
    for cc=1:cmax
        subplot(rmax, cmax, pIdx);   
   
        hold on;
        assinIdx = find(assignments == pIdx);
        A_binary_sub = A_binary(:, assinIdx);
        
        plot(A_binary_sub, 'lineWidth', 2);
        plot(mean(A_binary_sub'), 'k', 'lineWidth', 2);
        
        ylim([-0.1 1.1]);
        title(sprintf('%d-th, size:%d', pIdx, length(assinIdx)));
        hold off;
        pIdx = pIdx+1;
    end
end

end

