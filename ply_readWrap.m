function [ model ] = ply_readWrap ( Path )
 model = [];
 
 disp(sprintf('Reading... %s', Path));
 [tri, ~, data, ~] = ply_read(Path, 'tri');
 
 
 if isfield(data, 'vertex') == 1
    
     model.v = [data.vertex.x data.vertex.y data.vertex.z];
     model.f = tri';
     model.c = [data.vertex.red data.vertex.green data.vertex.blue];
     model.n = [data.vertex.nx data.vertex.ny data.vertex.nz];
     if isfield(data.vertex, 'radius') == 1
        model.radius = data.vertex.radius;
     end
 end
  
end