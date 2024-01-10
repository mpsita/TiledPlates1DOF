% Clear
close all, clearvars

%%
%  Input
nx=4; ny=8; % number of double tiles along the x and y directions


deg2rad = pi/180;
A = 120*deg2rad;
B = 120*deg2rad;
D = 2*pi - A - B;

C = 90*deg2rad;

a = 1; 
b = 0.5;    c = 0.5;
d = a;
e = c;


lb = - b * cos(A);
le = - e * cos(D);
la = a * cos(A+B+C);
lc = - c * cos(C);
lf = lc + la + le - a - lb;

hb = b * sin(A);
he = e * sin(D);
ha = a * sin(A+B+C);
hc = c * sin(C);
hf = hb+he+ha-hc;

f = (lf^2 + hf^2)^0.5;
F = pi-atan(hf/lf);
E = 4 * pi - A - B - C - D - F ;

abs(sin(A)*sin(E)/sin(B)/sin(F))


lf2 = f * cos(E-(pi-D));
hf2 = f * sin(E-(pi-D));
lb2 = b * cos(B-(pi-C));
hb2 = b * sin(B-(pi-C));


% shifts
transly = [-lf+lc-lb2+lc, hf+hc+hb2+hc];

% traslaz direz d1
tx1 = [a + f*cos(E-pi+C) + a + b*cos(pi-A), -f*sin(E-pi+C) + b*sin(pi-A)];

tx2 = [le + a + lc + la, - he + 0 + hc - ha ];
ty1 = [ lb-e*cos(D) , hb+e*sin(D) ];      ty2 = [ lf2-le , hf2+he ];

xynodes=[];

% nodes

% first row
% first double tile 
xynodes(1,:) = [0,0];
xynodes(2,:) = xynodes(1,:) + [a,0];
xynodes(3,:) = xynodes(2,:) + [f*cos(E-pi+C), -f*sin(E-pi+C)]; %%%%%%%%%%%%%%%%%%%%%%%
xynodes(4,:) = xynodes(3,:) + [a,0];
cn=4;
% next tiles
for ii=2:nx
    xynodes(cn+1:cn+4,:) = xynodes(cn-3:cn,:) + tx1;
    cn=cn+4;
end
%nr=[-lb, -hb];


% next rows
for jj=1:ny
    if mod(jj,2)==1
        % first tile
        xynodes(cn+1:cn+4,:) = [
            [e*cos(D) , e*sin(D)] ;
            [e*cos(D) , e*sin(D)]+[f*cos(E-pi+D),f*sin(E-pi+D)];
            [e*cos(D) , e*sin(D)]+[f*cos(E-pi+D),f*sin(E-pi+D)]+[a*cos(C-D),-a*sin(C-D)];
            [e*cos(D) , e*sin(D)]+[f*cos(E-pi+D),f*sin(E-pi+D)]+[a*cos(C-D),-a*sin(C-D)]+[b*sin(pi-A-(pi/2-C+D)),-b*cos(pi-A-(pi/2-C+D))];
            ];
        nr=xynodes(cn+1,:);
        cn=cn+4;
        % next tiles
        for ii=2:nx
            xynodes(cn+1:cn+4,:) = xynodes(cn-3:cn,:) + tx2;
            cn=cn+4;
        end
        % last side
        xynodes(cn+1,:)=xynodes(cn,:)+[la,-ha];
        xynodes(cn+2,:)=xynodes(cn+1,:)+[le,-he];
        cn=cn+2;
    else
        % first tile
        xynodes(cn+1:cn+4,:) = [
            nr+ty2;
            nr+ty2+[lb,hb];
            nr+ty2+[lb,hb]+[a,0];
            nr+ty2+[lb,hb]+[a,0]+[lf,-hf];
            ];
        nr=xynodes(cn+1,:);
        cn=cn+4;
        % next tiles
        for ii=2:nx
            xynodes(cn+1:cn+4,:) = xynodes(cn-3:cn,:) + tx1;
            cn=cn+4;
        end
        % last side
        xynodes(cn+1,:)=xynodes(cn,:)+[a,0];
        xynodes(cn+2,:)=xynodes(cn+1,:)+[le,-he];
        cn=cn+2;
    end
end


% tiles

tiles = [];

% first row
% first double tile
tiles(1,:)=[1 2 4*nx+4 4*nx+3 4*nx+2 4*nx+1];
tiles(2,:)=[2 3 4 5 4*nx+5 4*nx+4];
ct=2; % current tile
% next double tiles
for ii=1:nx-1
    tiles(ct+1,:)=tiles(ct-1,:)+4;
    tiles(ct+2,:)=tiles(ct,:)+4;
    ct=ct+2;
end
cs=8*nx+2;
tiles(ct,4)=cs;
% next rows
cn=4*nx;
for jj=2:ny
    for ii=1:nx
    tiles(ct+1,:)=[cn+2 cn+3 cn+nx*4+6 cn+nx*4+5 cn+nx*4+4 cn+nx*4+3];    
    tiles(ct+2,:)=[cn+3 cn+4 cn+5 cn+6 cn+nx*4+7 cn+nx*4+6];
    ct=ct+2;
    cn=cn+4;
    end
    tiles(ct,4)=cs+nx*4+2;
    cs=cs+4*nx+2;
    cn=cn+2;
end
% connections

connections=[];

% first row
% first one
connections(1,:)=[2 4*nx+4 1 2];
cc = 1;
ct = 1;
cn = 4;
% next ones
for ii = 2:nx
    connections(cc+1,:)=[cn+1 cn+4*nx+1 ct+1 ct+2];
    connections(cc+2,:)=[cn+2 cn+4*nx+4 ct+2 ct+3];
    cc=cc+2; ct=ct+2; cn=cn+4;
end
ct=ct+1;
% next rows
for jj = 2:ny
    cn=cn+1;
    % first one
    connections(cc+1,:)=[cn+1 cn+2 ct-1-(nx-1)*2 ct+1];
    connections(cc+2,:)=[cn+2 cn+4*nx+5 ct+1 ct+2];
    connections(cc+3,:)=[cn+2 cn+3 ct-1-(nx-1)*2 ct+2];
    connections(cc+4,:)=[cn+3 cn+4 ct-(nx-1)*2 ct+2];
    cn=cn+4; ct=ct+2; cc=cc+4;
    for ii = 2 : nx
        connections(cc+1,:)=[cn cn+1 ct-1-(nx-1)*2 ct];
    connections(cc+2,:)=[cn+1 cn+4*nx+2 ct ct+1];
    connections(cc+3,:)=[cn+1 cn+2 ct-1-(nx-1)*2 ct+1];
    connections(cc+4,:)=[cn+2 cn+4*nx+5 ct+1 ct+2]; 
    connections(cc+5,:)=[cn+2 cn+3 ct-1-(nx-1)*2 ct+2];
    connections(cc+6,:)=[cn+3 cn+4 ct-(nx-1)*2 ct+2];
    cc=cc+6; cn=cn+4; ct=ct+2;
    end
    cn=cn+1;
end


nodes=zeros(3,size(xynodes,1));
nodes(1:2,:)=xynodes';
edges=[]; stypes=[];
edges=connections(:,1:2);
stypes=ones(size(edges,1),1);
% firstpicture1_(nodes,edges,stypes);


%%
%  Creating files

fid = fopen('Nodes.dat','w');
fprintf(fid, '%f %f \n', xynodes.');
fclose(fid);

fileName = fopen('Tiles.dat','w');
fprintf(fileName, '%d %d %d %d %d %d\n', tiles.');
fclose(fileName);

fileName2 = fopen('Connections.dat','w');
fprintf(fileName2, '%d %d %d %d\n', connections.');
fclose(fileName2);
