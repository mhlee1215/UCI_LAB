function [ v, f, c ] = loadData3D( filename )

v = [];
f = [];
c = [];

    [pathstr,name,ext] = fileparts(filename); 

    if strcmpi(ext, '.obj') == 1
        obj = loadawobj(filename);
        v = obj.v';
        f = obj.f3';
    elseif strcmpi(ext, '.ply') == 1
        [tri, pts, data, comments] = ply_read(filename, 'tri');
        v = [data.vertex.x data.vertex.y data.vertex.z];
        f = tri';
        c = [data.vertex.red data.vertex.green data.vertex.blue];
    end


    

end

