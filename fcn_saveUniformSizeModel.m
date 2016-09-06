function [] = fcn_saveUniformSizeModel( Rv, Rc, Rn, dataRoot, fileName, density )
%FCN_SAVEUNIFORMSIZEMODEL Summary of this function goes here
%   Detailed explanation goes here
% fileName = 'merged_ref_big';
% density = 50;

if iscell(Rv)
    Rv2 = [];
    Rc2 = [];
    Rn2 = [];
    for i=1:length(Rv)
        Rv2 = [Rv2 Rv{i}];
        Rc2 = [Rc2 Rc{i}];
        if ~isempty(Rn)
            Rn2 = [Rn2 Rn{i}];
        end
    end
    
    if size(Rv2, 1) < size(Rv2, 2)
        Rv2 = Rv2';
        Rc2 = Rc2';
        Rn2 = Rn2';
    end
    
    fcn_saveUniformSizeModel( Rv2, Rc2, Rn2, dataRoot, fileName, density );
    
else
    if density == 0
        vv = Rv;
        cc = Rc;
        nn = Rn;
        writingPath = sprintf('%s/%s.ply', dataRoot, fileName);
    else
        [vv, cc, nn, ~] = fcn_uniformSampling(Rv', Rc', Rn', density);    
        writingPath = sprintf('%s/%s_d%d.ply', dataRoot, fileName, density);
    end

    d.vertex.x = vv(:,1);
    d.vertex.y = vv(:,2);
    d.vertex.z = vv(:,3);
    d.vertex.red = uint8(cc(:,1).*255);
    d.vertex.green = uint8(cc(:,2).*255);
    d.vertex.blue = uint8(cc(:,3).*255);
    if ~isempty(nn)
        d.vertex.nx = nn(:,1);
        d.vertex.ny = nn(:,2);
        d.vertex.nz = nn(:,3);
    end
    % d.vertex.radius = r;


    disp(sprintf('Writing Ply...%s', writingPath));
    % progressbar2(curStep/maxStep); curStep = curStep + 1;
    ply_write(d, writingPath, 'binary_little_endian');
end
% ply_write(d, '/home/mhlee/data_from_odroid/match/LAB_4_ref_tssss.ply', 'binary_little_endian');

% tic;
% ply_write(d, writingPath, 'ascii');
% toc;

end

