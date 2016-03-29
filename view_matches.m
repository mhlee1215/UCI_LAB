function [ output_args ] = view_matches( vertex11, color1, vertex22, color2, key11, key22, NN_Data, aa, validRatio )
%VIEW_MATCHES Summary of this function goes here
%   Detailed explanation goes here

if nargin < 9
    validRatio = 0.8;
end

showAll = true;

h=figure;
for i=1:size(key11, 2)/10
interval = 1;
clf;
hold on;

scatter3(vertex11(1:interval:end,1), vertex11(1:interval:end,2), ...
    vertex11(1:interval:end,3), 8, color1(1:interval:end, :), 'filled');

scatter3(vertex22(1:interval:end,1), vertex22(1:interval:end,2), ...
    vertex22(1:interval:end,3), 8, color2(1:interval:end, :), 'filled');

scatter3(key11(1,:)', key11(2,:)', key11(3,:)', '+');
scatter3(key22(1,:)', key22(2,:)', key22(3,:)', '*');

view(28, -20);
axis equal;


x = [key11(1, :)' key22(1, NN_Data(:,1))']';
y = [key11(2, :)' key22(2, NN_Data(:,1))']';
z = [key11(3, :)' key22(3, NN_Data(:,1))']';

% validIdx = find((aa(:,1) > 0).* (aa(:,1) ./ aa(:,2) < 0.7) .* (aa(:,1) < 0.03));

validIdx = find((aa(:,1) ./ aa(:,2) < validRatio));
range = 10*(i-1)+1:min(10*i, length(validIdx));
if showAll
    plot3(x(:,validIdx), y(:,validIdx), z(:,validIdx), 'lineWidth', 2);
    scatter3(x(1,validIdx), y(1,validIdx), z(1,validIdx), 200, 'filled');
    scatter3(x(2,validIdx), y(2,validIdx), z(2,validIdx), 200, 'filled');
    break;
else
    plot3(x(:,validIdx(range)), y(:,validIdx(range)), z(:,validIdx(range)), 'lineWidth', 2);
    scatter3(x(1,validIdx(range)), y(1,validIdx(range)), z(1,validIdx(range)), 200, 'filled');
    scatter3(x(2,validIdx(range)), y(2,validIdx(range)), z(2,validIdx(range)), 200, 'filled');
end
% saveas(h, sprintf('results/EF_match/matches_%d.png', i));
% saveas(h, sprintf('results/EF_match/matches_%d.fig', i));
pause();
end


end

