function [xyNodes, Tiles] = import_NT
% 
% Reads from file Nodes and Tiles
%
% Import Nodes data
f1 = fopen('Nodes.dat');
TxyNodes = fscanf(f1,'%g %g',[2 inf]);
fclose(f1);
xyNodes = TxyNodes';

%Import Tiles data
fid = fopen('Tiles.dat');
line1 = fgetl(fid);
res=line1;
while ischar(line1)
line1 = fgetl(fid);
res =char(res,line1);
end
fclose(fid);
for k=1:size(res,1)-1
  Tiles{k}=str2num(res(k,:));
end
