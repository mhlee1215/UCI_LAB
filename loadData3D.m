function [ v, f, c, g ] = loadData3D( filename )

v = [];
f = [];
c = [];
g = [];

    [pathstr,name,ext] = fileparts(filename); 

    if strcmpi(ext, '.obj') == 1
        obj = loadawobj(filename);
        v = obj.v';
        f = obj.f3';
        g = obj.g3;
    elseif strcmpi(ext, '.ply') == 1
        [tri, pts, data, comments] = ply_read(filename, 'tri');
        v = [data.vertex.x data.vertex.y data.vertex.z];
        f = tri';
        c = [data.vertex.red data.vertex.green data.vertex.blue];
    end


    

end

