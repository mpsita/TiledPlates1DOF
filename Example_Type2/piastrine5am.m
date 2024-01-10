% 
% MATLAB code: piece-wise rigid plates
%

%%
% clear memory
clc
close all
clearvars

%%
% Input Nodes and Tiles data from file
[xyNodes, Tiles] = import_NT;
%%
% Input Nodes and Tiles data manually
% 
% Example 1
%xyNodes = [-1 -1;0 0; 0 1 ;1 -1; 1 -2; 2 0];                      % Coordinates of the nodes
%Tiles = {[1 2 3]; [1 2 4 5]; [3 2 4 6]; [5 4 6]};                 % Tiles defined by nodes
% Example 2
%xyNodes = [0 0;1 0; 2 0;3 0; 1 -1; 2 -1; 1 1; 2 1];                   % Coordinates of the nodes
%Tiles = {[1 2 7]; [1 2 3 6 5]; [2 3 8 7]; [3 6 4 8]};                 % Tiles defined by nodes
% Example 3  two-degrees of freedom
%xyNodes = [-1 -1;0 -1; 1 -1;-1 0; 0 0; 1 0; -1 1; 0 1; 1 1];                   % Coordinates of the nodes
%Tiles = {[1 2 5 4]; [2 3 6 5]; [4 5 8 7]; [5 6 9 8]};                 % Tiles defined by nodes
%%
% Size of the input
Nnodes = size(xyNodes,1);                      % Number of nodes
Ntiles = length(Tiles);                        % Number of tiles
Nposconnections=Nnodes*(Nnodes+1)/2-Nnodes;    % Number of possible connections
SizeTile = cellfun('length',Tiles);            % Vector whose i component contains dimension of tile i
%%
% Drawing tiles
FV.vertices=xyNodes;
figure(1);
for i = 1:Ntiles
    FV.faces=Tiles{i};
    XX = xyNodes(FV.faces,:);
    CG=sum(XX)/6;
    patch(FV,'facecolor','y', 'FaceAlpha', 0.2);
    text(CG(1),CG(2),num2str(i));
    
    for j=1:6
      text(XX(j,1),XX(j,2),num2str(FV.faces(j)));  
    end
end



%%
fprintf('------------------------------------\n')
fprintf('Number of nodes = %d \n', Nnodes)
fprintf('Number of tiles = %d \n', Ntiles)

if 1
% Evaluate Connections
tic
Connections=eval_connections(Nposconnections,Tiles);
toc
%Connections=eval_connections_old(Nposconnections,Tiles)
else
% Import Connections data from file
fc1 = fopen('Connections.dat');
Connections = fscanf(fc1,'%d %d %d %d',[4 inf]);
fclose(fc1);
Connections = Connections.';
end

%%%%
Nconstraints=length (Connections(:,1));
fprintf('Number of constraints = %d \n', Nconstraints)

%%
% Kinematical matrix
%
K=zeros(Nconstraints,2*Ntiles);
for i=1:Nconstraints
    xdiff= xyNodes(Connections(i,1),1)-xyNodes(Connections(i,2),1);
    ydiff= xyNodes(Connections(i,1),2)-xyNodes(Connections(i,2),2);
    K(i,(Connections(i,3)-1)*2+1)=xdiff;
    K(i,(Connections(i,3)-1)*2+2)=ydiff;
    K(i,(Connections(i,4)-1)*2+1)=-xdiff;
    K(i,(Connections(i,4)-1)*2+2)=-ydiff;
end
rk=rank(K,10^-3);
fprintf('The kinematic matrix has dimensions %d times %d and its rank is equal to %d\n', Nconstraints, 2*Ntiles,rk)
%%
% Evaluation of infinitesimal rotations
%
fprintf('------------------------------------\n')
fprintf('To avoid infinitesimal rigid rotations set the degrees of freedom of one tile equal to zero\n');
fprintf('Choose tile (between 1 and %d)', Ntiles);
Tzero =input(': ');
Kred=K;
Kred(:,2*Tzero-1:2*Tzero)=[];   %reduced kinematic matrix (without zero rotations)
[a,b,c]=svd(Kred);
Zred=c(:,rk+1:end);                % null space of the reduced matrix
svals=diag(b);
fprintf('last singular values:\n')
fprintf(num2str(svals(rk-2:end)'))
fprintf('\n')
if rk<min(size(Kred))
    svr=svals(rk+1)/svals(rk);
    fprintf('Singular values ratio: %d \n', svr)
end
dof=2*Ntiles-rk-2
doff=size(Kred,2)-rk;
Zzero=zeros(2,dof);
Z = [Zred(1:(Tzero-1)*2,:);Zzero;Zred((Tzero-1)*2+1:end,:)];    % null space of the kinematic matrix
if dof==1
    fprintf('There is just 1 degree of freedom.\n')
else
    fprintf('There are %d degrees of freedom.\n', dof)
end

%%
% Evaluation of infinitesimal displacements
%
fprintf('------------------------------------\n')
fprintf('To avoid infinitesimal translations set the displacement of a node (between 1 and %d) equal to zero\n', Nnodes);
Nzero =input('Choose node: ');
u = zeros(Nnodes, dof);

for k=1:dof
Tcheck=zeros(Ntiles,1);
Ncheck=zeros(Nnodes,1);
Ncheck(Nzero)=1;
while sum(Ncheck)<Nnodes
  for r=1:Nnodes
    if Ncheck(r)==1
      for j=1:Ntiles
         if ismember(r,Tiles{j}) && Tcheck(j)==0 
            Tcheck(j)=1; 
            for i=1:SizeTile(j) 
                 if Ncheck(Tiles{j}(i))== 0
                     u(Tiles{j}(i),k)=u(r,k)+Z(2*j-1,k)*(xyNodes(Tiles{j}(i),1)-xyNodes(r,1))+ Z(2*j,k)*(xyNodes(Tiles{j}(i),2)-xyNodes(r,2));
                     Ncheck(Tiles{j}(i))=1;
                 end    
            end
         end
      end
    end
  end
end
end  % k cycle
fprintf('\n');
for rr = 1 : Nnodes
    fmt=['u%d =' repmat(' %d \t',1,numel(u(rr,:)))];
    fprintf(fmt,rr, u(rr,:));
    fprintf('\n');
end
%%
%  Plot displacement
%

fprintf('------------------------------------\n')
fprintf('To plot, set the values of the %d degrees of freedom. \n', dof);
Amp=zeros(dof,1);
for k=1:dof
    fprintf('Amplificator factor %d', k);
    Amp(k) =input(': ');
end
pu=zeros(Nnodes,1);
for k=1:dof
    pu=pu+Amp(k)*u(:,k);
end

xyzNodes=[xyNodes,pu];
FV.vertices=xyzNodes;
figure(2);
for i = 1:Ntiles
    faces=Tiles{i};
    patch('Faces', faces, 'Vertices', xyzNodes, 'FaceColor', "yellow"); %"#EDB120");
end
% for i = 1:Ntiles
%     faces=Tiles{i};
%     patch('Faces', faces, 'Vertices', xyNodes, 'FaceColor', "none");
% end
% Axes settings
xlabel('x'); ylabel('y'); zlabel('z');
axis vis3d equal;
view([70, 15]);
camlight;
set(gca,'XColor', 'none','YColor','none', 'ZColor', 'none')
set(gca, 'color', 'white');
grid off





if dof>1
fprintf('------------------------------------\n')
fprintf('Do you want to change the amplificator factors? \n', dof);
Repeat =input('Type 1 if you want to change them, 0 otherwise: ');
end

while dof>1 && Repeat==1
fprintf('------------------------------------\n')
fprintf('To plot, set the values of the %d degrees of freedom. \n', dof);
Amp=zeros(dof,1);
for k=1:dof
    fprintf('Amplificator factor %d', k);
    Amp(k) =input(': ');
end
pu=zeros(Nnodes,1);
for k=1:dof
    pu=pu+Amp(k)*u(:,k);
end

xyzNodes=[xyNodes,pu];
FV.vertices=xyzNodes;
figure(2);
for i = 1:Ntiles
    faces=Tiles{i};
    patch('Faces', faces, 'Vertices', xyzNodes, 'FaceColor', 'y');
end
% Axes settings
xlabel('x'); ylabel('y'); zlabel('z');
axis vis3d equal;
view([30, 30]);
camlight;
grid on;

fprintf('------------------------------------\n')
fprintf('Do you want to change the amplificator factors? \n', dof);
Repeat =input('Type 1 if you want to change them, 0 otherwise: ');
end  %while