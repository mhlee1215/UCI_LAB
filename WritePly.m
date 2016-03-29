%Develop by Renoald
%University Teknologi Malaysia
%For academic purpose
%Email:renoald@live.com
%This function write ply file
%mention:This function only read ply ascii file 
%since for acedemic purpose , this file not support color and
%other option for ply format
function WritePly(p, f, file)
fid=fopen(file,'w');
% f=dlmread('face');
% p=dlmread('point');
fprintf(fid,'ply\n');
fprintf(fid,'format ascii 1.0\n');
fprintf(fid,'comment VCGLB generated\n');
fprintf(fid,'element vertex %d\n',length(p));
fprintf(fid,'property float x\n');
fprintf(fid,'property float y\n');
fprintf(fid,'property float z\n');
fprintf(fid,'element face %d\n',length(f));
fprintf(fid,'property list uchar int vertex_indices\n');
fprintf(fid,'end_header\n');
fid1=fopen('point');
while 1
    tline = fgetl(fid1);
    if ~ischar(tline),   break,   end
    disp(tline)
    fprintf(fid,'%s\n',tline);
end
fclose(fid1)
fid2=fopen('face');

while 1
    tline1 = fgetl(fid2);
    if ~ischar(tline1),   break,   end
    disp(tline1)
    fprintf(fid,'%s\n',tline1);
end
fclose(fid2);
fclose(fid);
