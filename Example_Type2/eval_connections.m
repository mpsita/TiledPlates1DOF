function [connections] = eval_connections(dummy,Tiles)
%	Summary of this function goes here
%   Detailed explanation goes here
edgetile=[];
Ntiles=size(Tiles,2);
% list edges and corresponding tiles
for ii=1:Ntiles
    for jj=1:size(Tiles{ii},2)-1
        edgetile=[edgetile; [Tiles{ii}(jj:jj+1) ii]];
    end
edgetile=[edgetile; [Tiles{ii}(end) Tiles{ii}(1) ii]];
edges=edgetile(:,1:2);

end

edges=sort(edges')';

% eliminate duplicates, track tiles
[C,IA,IC] = unique(edges,'rows','first');
[C2,IA2,IC2] = unique(edges,'rows','last');


connections=[C edgetile(IA,3) edgetile(IA2,3) ];

% eliminate self connections (boundary edges)
elimcon=[connections(:,3)==connections(:,4)];
connections(elimcon,:)=[];


