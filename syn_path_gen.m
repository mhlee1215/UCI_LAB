function [ poseSet ] = syn_path_gen( p, dt, vertices, type, pc )
%SYN_PATH_GEN Summary of this function goes here
%  in : N x 3 Mat



% dt = 0.01;
curve = p';
% interpote a 3d curve using spline 
% path 3*x 
% newPath 3*x

x = curve(1, :); 
y = curve(2, :); 
z = curve(3, :);

t = cumsum([0;sqrt(diff(x(:)).^2 + diff(y(:)).^2 + diff(z(:)).^2)]); 
sx = spline(t,x); 
sy = spline(t,y); 
sz = spline(t,z);

tt = t(1):dt:t(end); 
xp = ppval(sx, tt); 
yp = ppval(sy, tt); 
zp = ppval(sz, tt);

newCurve = [xp; yp; zp]; 
pp = newCurve';

% pc = mean(pp);
% pc = [0 1.6 0];

if ~isempty(vertices)
    figure; hold on;
    scatter3(x, y, z, 100, 'o', 'filled');
    scatter3(x(1), y(1), z(1), 50, 'r', 'filled');
    scatter3(x(end), y(end), z(end), 50, 'g', 'filled');
    scatter3(pp(:,1), pp(:,2), pp(:,3), 10, 'b', 'filled');
    scatter3(vertices(1:10:end,1),vertices(1:10:end,2), vertices(1:10:end,3), 1.5, 'filled');
    scatter3(pc(1), pc(2), pc(3), 100, 'r', 'filled');
    axis equal;
    grid;
end

poseSet = {};

% pc = mean(pp);
for i=1:size(pp, 1)-1
    if type == 1
        p1 = pp(i,:);
        p2 = pp(i+1,:);
      
        r = vrrotvec([0 0 -1], -(p2 - p1));
        m = vrrotvec2mat(r);
        
        m = m.*[1 1 1 ; -1 -1 -1; 1 1 1];
        zz = eye(3);

        m = inv(m)*zz;
        poseSet{end+1}.R = m;
        poseSet{end}.t = -m*p1';
    
    elseif type == 2
        po = pp(i,:);
        p1 = po;
        
        p1 = p1 - pc;
%         p2 = p1.*2;
      
         base = [0 0 1];
%        
%         r = RotMat([p1(1) p1(2) p1(3)], base);
%         m2 = r';
        
        r1 = RotMat(base, [p1(1) 0 p1(3)]);
%         r2 = RotMat(r1*p1', base);
        
%         r1*p1'
        
        m = r1;%r2';
        
        r2 = RotMat(m*[0 0 1]', [p1(1) p1(2) p1(3)]);
        m2 = r2*m;
        
        
        %Perspective correction - constratin y for upright.
        pu = p1 ./ sqrt(sum(p1.^2));
        
        if p1(2) >= 0
            t = angle3(pu, [0 1 0]);
            ux = 1/cos(t);

            up = [0 ux 0]-pu;
            up = up ./ sum(sqrt(up.^2));
        else
            t = angle3(pu, [0 -1 0]);
            ux = 1/cos(t);

            up = [0 ux 0]+pu;
            up = up ./ sum(sqrt(up.^2));
        end
        
        oldUp1 = (m2*[0 1 1]'-m2*[0 0 1]')';
%         oldUp2 = (m2*[0 1 1]'-m2*[0 0 1]')';
%         

%         r3 = RotMat(up, oldUp1)    
%         r4 = RotMat(oldUp2, oldUp1);       
        
        
%         m2 = r3*m2;
        
%         r3 = RotMat(m2*[0 1 0]', [0 1 0]);
%         m2 = r3*m2;
%         m = r2';
        
%         r2 = vrrotvec([p1(1) 0 p1(3)], [p1(1) p1(2) p1(3)]);
%         r2 = vrrotvec([0 0 1], [p1(1) p1(2) p1(3)]);
%         m2 = vrrotvec2mat(r2);
%         
%         m = m * m2;
        
%         m = m.*[1 1 1 ; -1 -1 -1; 1 1 1];
        
%         m = [1 0 0 ; 0 -1 0 ; 0 0 1];%rotationmat3D(pi/2, p1);
        zz = eye(3);
%         zz(2) = -1;
        
        m3 = inv(m2)*zz;

        
%         oldUp = (m3'*[0 1 1]'-m3'*[0 0 1]')';
        
%         m3 = eye(3);
        poseSet{end+1}.R = m3;
        poseSet{end}.t = -m3*po';
        poseSet{end}.up = up;
        poseSet{end}.oldUp = oldUp1;
        
    end
    
end

if ~isempty(vertices)
    imW = 1;
    imH = 1;
    for i=1:1:0%length(poseSet)-10
        cameraMat = poseSet{i};
    %    pM(:,end)
    %     scatter3(cameraMat.t(1), cameraMat.t(2), cameraMat.t(3), 0.1);
    %     P = projectImg2World(pM, [0 0]', 0.1);

        cameraMat.f = 1;
        cameraMat.imSize = [0 0];%[imH imW];

        up = cameraMat.up;
        oldUp = cameraMat.oldUp;
    %     cameraMat.R
        Cp = -cameraMat.R'*cameraMat.t;

        depth = 0.5;
        C1 = projectImg2World(cameraMat, [-imH -imW]', depth);
        C2 = projectImg2World(cameraMat, [-imH imW]', depth);
        C3 = projectImg2World(cameraMat, [imH imW]', depth);
        C4 = projectImg2World(cameraMat, [imH -imW]', depth);
        Cc = projectImg2World(cameraMat, [0 0]', depth);
        Cup = projectImg2World(cameraMat, [0 imH]', depth);
        PScene = [C1 ; C2 ; C3 ; C4 ; C1];

    %     scatter3(Cp(1), Cp(2), Cp(3), 50, 'm', 'filled');
        for i=1:size(PScene, 1)-1
            line([Cp(1) PScene(i,1)]', [Cp(2) PScene(i,2)]', [Cp(3) PScene(i,3)]', 'linewidth', 1);    
            line([PScene(i,1) PScene(i+1,1)]', [PScene(i,2) PScene(i+1,2)]', [PScene(i,3) PScene(i+1,3)]', 'linewidth', 1);    
            if i == 1
                color = 'g';
            elseif i == 2
                color = 'b';
            elseif i == 3
                color = 'r';
            elseif i == 4
                color = 'y';
            else color = 'k';

            end
            scatter3(PScene(i,1), PScene(i,2), PScene(i,3), 40, color, 'filled');
        end
        line([Cc(1) Cp(1)]', [Cc(2) Cp(2)]', [Cc(3) Cp(3)]', 'linewidth', 3);    
%         line([Cc(1) Cup(1)]', [Cc(2) Cup(2)]', [Cc(3) Cup(3)]', 'linewidth', 6, 'color', 'r');    
%         line([Cc(1) Cc(1)+up(1)]', [Cc(2) Cc(2)+up(2)]', [Cc(3) Cc(3)+up(3)]', 'linewidth', 5, 'color', 'g');    
%         line([Cc(1) Cc(1)+oldUp(1)]', [Cc(2) Cc(2)+oldUp(2)]', [Cc(3) Cc(3)+oldUp(3)]', 'linewidth', 3, 'color', 'k');    
        line([pc(1) Cp(1)]', [pc(2) Cp(2)]', [pc(3) Cp(3)]', 'linewidth', 1);    

    end

end




end

