clear all 
close all

x0      = [0, 0, 0];    %location of lower left edge of hypercube
expand  = 1;            %expand hypercube
Nrefine = 0;            %number of mesh-refinements
Nq      = 5;            %number of querypoints

%% create mesh
mesh = ofem_v2.Geometry;           
mesh.hypercube(x0,expand);


for n = 1:Nrefine
    
    mesh.uniform_refine;
    
end

mesh.create_edges;
co = double(reshape(permute(mesh.co,[3,1,2]),[],size(mesh.co,1)));

%% create querypoints and pointlocation
xq = [expand*rand(Nq,1), expand*rand(Nq,1), expand*rand(Nq,1)];
[idx, tr, bary] = ofem_v2.tools.pointLocation(mesh,xq,[]);
        
%% plot solution
for i = 1:numel(idx)
    
    figure(i)
    tetramesh(mesh.el,co,'FaceColor',[0 0 1],'FaceAlpha',0.1);
    hold on 
    plot3(xq(i,1),xq(i,2),xq(i,3),'x','MarkerSize',9,'LineWidth',2)
    tetramesh(mesh.el(idx(i),:),co,'FaceColor',[1 0 0],'FaceAlpha', 0.2);
    hold off
 
end