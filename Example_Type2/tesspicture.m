function tesspicture(xynodes,tiles,connections)
% this function draw a picture of the structure after data generation

nodes=zeros(3,size(xynodes,1));
nodes(1:2,:)=xynodes';
edges=connections(:,1:2);

Nedges=(size(edges,1));
Nnodes=(size(nodes,2));

l=length0(nodes,edges);

f=figure;

grid off;

% handles of nodes
q=zeros(Nnodes,1);
% handles of edges
p=zeros(Nedges,1);

%%%%%%%%%%%%%%%%%%%
%%%%%% nodes %%%%%%
%%%%%%%%%%%%%%%%%%%

%for k=1:Nnodes
%         q(k)=draw_node(nodes(:,k));
%end

nt=zeros(Nnodes,1);
dt=0;
%dt=max(nodes(:))/100; % text distance from a node
for k=1:Nnodes
         nt(k)=text_node(k,dt,nodes(:,k));
    %              set(nt(k),'FontSize',11,'FontWeight','normal');
end

for k=1:size(tiles,2)
    center=sum(nodes(:,tiles{k})')'/size(tiles{k},2);
         nt(k)=text_node(k,dt,center);
         set(nt(k),'FontSize',11,'FontWeight','normal');
end

%%%%%%%%%%%%%%%%%%%
%%%%%% edges %%%%%%
%%%%%%%%%%%%%%%%%%%

for k=1:Nedges
         p(k)=draw_line(l(k),dt,nodes(:,edges(k,1)),nodes(:,edges(k,2)));
end


for k=1:Nedges
        set(p(k),'LineWidth',1);
end

% et=zeros(Nedges,1);
% for k=1:Nedges
%          et(k)=text_line(k,dt,etypes(k),nodes(:,edges(k,1)),nodes(:,edges(k,2)));
% end

axis equal
axis on, grid off
set(gcf,'color','white')

br=1.1;
gap=0.1;
x1=br*min(nodes(1,:))-0.01; x2=br*max(nodes(1,:))+0.01;
x1=x1-gap*(x2-x1);       x2=x2+gap*(x2-x1);
y1=br*min(nodes(2,:))-0.01; y2=br*max(nodes(2,:))+0.01;
y1=y1-gap*(y2-y1);       y2=y2+gap*(y2-y1);
z1=br*min(nodes(3,:))-0.01; z2=br*max(nodes(3,:))+0.01;
z1=z1-gap*(z2-z1);       z2=z2+gap*(z2-z1);
set(gca,'Zlim',[z1 z2],'Xlim',[x1 x2],'Ylim',[y1 y2]);

X=text('string','X'); Y=text('string','Y'); Z=text('string','Z');
set(gca,'xlabel',X,'ylabel',Y,'zlabel',Z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function n=draw_node(coords)

n=line(coords(1),coords(2),coords(3),...
    'Marker','o',...
    'MarkerSize',12,...
    'MarkerEdgeColor','red',...
	'MarkerFaceColor','yellow');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function nt=text_node(k,dt,coords)
nt=text(coords(1),coords(2),coords(3),int2str(k),...
    'HorizontalAlignment','Center',...
    'FontSize',14,'FontWeight','Bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function p=draw_line(label,dt,coords1,coords2)
d=0.05;
p1=(1-d)*coords1+d*coords2;
p2=d*coords1+(1-d)*coords2;

a=[p1(1) p2(1)];
b=[p1(2) p2(2)];
c=[p1(3) p2(3)];


p=line(a,b,c,'color','blue','LineWidth',1);

uicm=uicontextmenu;
set(p,'uicontextmenu',uicm);
uimenu(uicm,'label',['L=',num2str(label)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function nt=text_line(k,dt,e,coords1,coords2)
coords=0.5*(coords1+coords2);

%line(coords(1),coords(2),coords(3),...
%    'Marker','Square',...
%    'MarkerSize',14,...
%    'MarkerEdgeColor','Blue',...
%	'MarkerFaceColor','yellow');
    
nt=text(coords(1),coords(2),coords(3),int2str(k),...
    'FontSize',9,...
    'Color','Black',...
    'HorizontalAlignment','Center',...
    'Visible','off'...
);
