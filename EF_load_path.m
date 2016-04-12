function [ poseSet ] = EF_load_path(filePath, vertices)

poseFilePath = filePath;%'syn_path/out.txt';
poseFileID = fopen(poseFilePath, 'r');
poseMat = fscanf(poseFileID, '%f');
poseSet = {};
imW = 3;
imH = 3;
for pii = 1:size(poseMat, 1)/19
    idx = (pii-1)*19+1;
    pM = reshape(poseMat(idx+3:idx+18)', 4, 4)';
    cameraMat.R = pM(1:3, 1:3);
    cameraMat.t = pM(1:3, end);
    cameraMat.imSize = [0 0];%[imH imW];
    cameraMat.f = 1;
    poseSet{end+1} = cameraMat;
end



if exist('vertices', 'var')
    h=figure; hold on;
    scatter3(vertices(1:10:end,1),vertices(1:10:end,2), vertices(1:10:end,3), 0.5, 'filled');
    
    poseMM = [];
    for i=1:length(poseSet)
        cameraMat = poseSet{i};
        camPos = -cameraMat.R'*cameraMat.t;
        poseMM = [poseMM ; camPos(1) camPos(2) camPos(3)];

    end
    scatter3(poseMM(:,1), poseMM(:,2), poseMM(:,3), 0.5);
    axis equal;
end

for i=1:100:length(poseSet)
    cameraMat = poseSet{i};
%    pM(:,end)
%     scatter3(cameraMat.t(1), cameraMat.t(2), cameraMat.t(3), 0.1);
%     P = projectImg2World(pM, [0 0]', 0.1);
    
    
    cameraMat.R
    Cp = -cameraMat.R'*cameraMat.t;
    depth = 0.5;
    C1 = projectImg2World(cameraMat, [-imH -imW]', depth);
    C2 = projectImg2World(cameraMat, [-imH imW]', depth);
    C3 = projectImg2World(cameraMat, [imH imW]', depth);
    C4 = projectImg2World(cameraMat, [imH -imW]', depth);
    PScene = [C1 ; C2 ; C3 ; C4 ; C1];

    for i=1:size(PScene, 1)-1
        line([Cp(1) PScene(i,1)]', [Cp(2) PScene(i,2)]', [Cp(3) PScene(i,3)]', 'linewidth', 1);    
        line([PScene(i,1) PScene(i+1,1)]', [PScene(i,2) PScene(i+1,2)]', [PScene(i,3) PScene(i+1,3)]', 'linewidth', 1);    
    end
    
end
axis equal;
saveas(h, sprintf('%s/trajectory_1.fig', resultPath));
close(h);
end
