clear all 
close all

x0      = [0, 0];   %location of lower left edge of hypercube
expand  = 1;        %expand hypercube
Nrefine = 2;        %number of mesh-refinements
Nq      = 5;        %number of querypoints

%% create mesh
mesh = ofem_v2.Geometry;
mesh.hypercube(x0,expand);


for n = 1:Nrefine
    
    mesh.uniform_refine;
    
end

mesh.create_edges;
co = double(reshape(permute(mesh.co,[3,1,2]),[],size(mesh.co,1)));

%% create querypoints and pointlocation
xq = [expand*rand(Nq,1), expand*rand(Nq,1)];
[idx, tr, bary] = ofem_v2.tools.pointLocation(mesh,xq,[]);
        
%% plot solution
for i = 1:numel(idx)

    figure(i)
    trimesh(mesh.el,co(:,1),co(:,2),'Color','k');
    hold on 
    plot(xq(i,1),xq(i,2),'rx','MarkerSize',9,'LineWidth',1);
    trimesh(mesh.el(idx(i),:),co(:,1),co(:,2),'Color','g');
    hold off

end