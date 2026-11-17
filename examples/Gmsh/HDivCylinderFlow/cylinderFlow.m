%% Potential flow past a cylinder in a channel with H(div) elements
%
% Mixed formulation of potential (Darcy) flow
%
%       u + grad p = f      in the channel
%            div u = 0      in the channel
%              u.n = 0      on the walls and the cylinder
%                p = 0      at the inlet and the outlet
%
% with the constant driving force f = (G,0). This is the same as
% prescribing the pressure drop G*L between inlet and outlet: the physical
% pressure is P = p - G*x. The velocity is irrotational and divergence
% free, i.e. potential flow around the cylinder.
%
% Discretisation:
%   velocity  u in HDivElement(2,order): BDM_order (Raviart-Thomas for 0)
%   pressure  p discontinuous P_{order-1} (P_0 for order 0)
% Since div(BDM_k) = P_{k-1}, the discrete velocity is exactly divergence
% free, pointwise and on every element.
%
% Weak form: find (u,p) with
%   (u,v) - (p,div v) = (f,v)   for all v in V_h with v.n = 0 on walls
%           (div u,q) = 0       for all q in Q_h
% The pressure boundary condition p = 0 is natural (no boundary term).

close all
clear

order = 2;      % polynomial order of the velocity
h     = 1;      % mesh refinement (mesh size ~ 0.04/h)
G     = 1;      % driving pressure gradient

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
file = fullfile(here,'channel');

%% Mesh
[status,msg] = system(sprintf('gmsh "%s.geo" -setnumber h %g -2 -v 0 -o "%s.msh"',file,h,file));
if status ~= 0
    error('gmsh failed: %s',msg);
end
mesh = ofem_v2.Geometry;
mesh.load_from_msh(file);
mesh.reorderAC;
mesh.create_edges;
mesh.jacobiandata;

fluid = ofem_v2.materials.Material;
fluid.mass = 1;                 % inverse permeability
mesh.setMaterial('Domain',fluid);

%% Velocity space and boundary conditions
fe = ofem_v2.elements.loadFE(sprintf('HDiv_2D_Order_%d',order));

dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(fe);

walls = ofem_v2.boundary.NormalTrace('Wall',mesh);
cylinder = ofem_v2.boundary.NormalTrace('Cylinder',mesh);

phys = ofem_v2.Physical_Problem(fe,mesh,0,1,0);
phys.setParaM('mass');
phys.setBoundaryCondition(walls);
phys.setBoundaryCondition(cylinder);
phys.attachDOFHandler(dofs);
dofs.generateDOFs;
phys.assemble;

M = phys.M;
F = fe.volumeForce(phys,[G;0],1:mesh.Nint);

%% Pressure space and divergence matrix B(i,j) = (div v_j, q_i)
kp = max(order-1,0);
B = divergenceMatrix(fe,phys,kp);
Np = size(B,1);

%% Saddle point system
free = dofs.freeDOFs;
Nf = numel(free);
A = [M(free,free), -B(:,free)'; -B(:,free), sparse(Np,Np)];
x = A\[F(free); zeros(Np,1)];

u = zeros(dofs.Nd,1);
u(free) = x(1:Nf);
p = x(Nf+1:end);

%% Evaluate velocity and pressure
[w,l] = ofem_v2.tools.gaussSimplex(2,2*order+2);
[X,U,divU] = fe.evaluate(phys,u,l);
P = evalPressure(p,kp,l,mesh.Nint) - G*X(1,:,:);   % physical pressure

detA = abs(reshape(mesh.detD,1,[]));                 % 2*|element area|
integrate = @(f) sum(sum(reshape(f,mesh.Nint,[]).*w(:)',2).*detA(:));

Lc = 2.2; Hc = 0.41;
Q = integrate(U(1,:,:))/Lc;                          % mean flow rate
% speed half way between the top of the cylinder and the wall. (The
% maximum speed sits on the cylinder, which the mesh approximates by a
% polygon; its corners make the discrete maximum grow with the order.)
xs = [0.2;0.33];
us = pointValue(fe,phys,u,xs);

fprintf('H(div) order %d (velocity DOFs %d, pressure DOFs %d)\n',order,Nf,Np);
fprintf('  max |div u_h| at quadrature points : %.2e\n',max(abs(divU(:))));
fprintf('  flow rate Q                         : %.6f (without cylinder %.6f)\n',Q,G*Hc);
fprintf('  speed at (%.2f,%.2f) / inflow speed : %.5f\n',xs,norm(us)/(Q/Hc));
fprintf('  pressure drop inlet -> outlet       : %.4f (prescribed %.4f)\n', ...
    max(P(:))-min(P(:)),G*Lc);

%% Plot speed, streamlines and pressure
co = reshape(mesh.co,2,[])';
tri = mesh.el(:,1:3);
Ue = squeeze(sum(U.*reshape(w,1,1,[]),3))*2;         % element mean velocity
Pe = squeeze(sum(P.*reshape(w,1,1,[]),3))*2;         % element mean pressure
cen = squeeze(mean(reshape(mesh.co(:,:,tri'),2,3,[]),2));

figure('Name','Potential flow past a cylinder','Position',[100 100 1100 600]);
subplot(2,1,1)
patch('Faces',tri,'Vertices',co,'FaceVertexCData',vecnorm(Ue,2,1)', ...
    'FaceColor','flat','EdgeColor','none');
hold on
[gx,gy] = meshgrid(linspace(0,Lc,441),linspace(0,Hc,83));
gu = scatteredInterpolant(cen(1,:)',cen(2,:)',Ue(1,:)','linear','nearest');
gv = scatteredInterpolant(cen(1,:)',cen(2,:)',Ue(2,:)','linear','nearest');
GU = gu(gx,gy); GV = gv(gx,gy);
inCyl = (gx-0.2).^2+(gy-0.2).^2 < 0.05^2;
GU(inCyl) = NaN; GV(inCyl) = NaN;
sy = linspace(0.012,Hc-0.012,15);
hs = streamline(gx,gy,GU,GV,0.005*ones(size(sy)),sy);
set(hs,'Color','k','LineWidth',0.75);
axis equal tight; colorbar; title('|u_h| and streamlines')
subplot(2,1,2)
patch('Faces',tri,'Vertices',co,'FaceVertexCData',Pe', ...
    'FaceColor','flat','EdgeColor','none');
axis equal tight; colorbar; title('pressure P = p_h - G x')

%% Local functions

function B = divergenceMatrix(fe,phys,k)
    % B(i,j) = int q_i div(v_j) with q_i the monomials of degree <= k in
    % the reference coordinates of each element (discontinuous P_k)
    g = phys.geometry;
    Ne = g.Nint;
    E = exponents(k);
    nq = size(E,1);
    [w,l] = ofem_v2.tools.gaussSimplex(2,max(k+max(fe.degree-1,0),1));
    Ns = fe.DOFsPerElement;
    Bl = zeros(nq,Ns,Ne);
    for q = 1:numel(w)
        dphi = cat(3,fe.divN{1}(l(1,q),l(2,q)),fe.divN{2}(l(1,q),l(2,q)));
        dphi = dphi(:,:,g.refTet);
        mq = prod(l(:,q)'.^E,2);
        Bl = Bl + w(q)*mq.*dphi;
    end
    % div v = div v_ref/detD and dx = |detD| dx_ref
    Bl = Bl.*sign(g.detD);
    rows = reshape(1:nq*Ne,nq,Ne);
    I = repmat(reshape(rows,nq,1,Ne),1,Ns,1);
    J = repmat(reshape(phys.DOFs.el2DOF',1,Ns,Ne),nq,1,1);
    B = sparse(I(:),J(:),Bl(:),nq*Ne,phys.DOFs.Nd);
end

function v = pointValue(fe,phys,u,x)
    % velocity at the point x (2 x 1)
    g = phys.geometry;
    x1 = g.co(:,:,g.el(:,1));
    lam = pagemtimes(pageinv(g.Dk),x-x1);            % reference coordinates
    lam = reshape(lam,2,[]);
    e = find(all(lam >= -1e-12,1) & sum(lam,1) <= 1+1e-12,1);
    [~,U] = fe.evaluate(phys,u,lam(:,e));
    v = U(:,e);
end

function P = evalPressure(p,k,l,Ne)
    % values of the discontinuous pressure at the reference points l
    E = exponents(k);
    nq = size(E,1);
    pe = reshape(p,nq,Ne);
    P = zeros(1,Ne,size(l,2));
    for q = 1:size(l,2)
        P(1,:,q) = prod(l(:,q)'.^E,2)'*pe;
    end
end

function E = exponents(k)
    E = zeros(0,2);
    for a = 0:k
        for b = 0:k-a
            E(end+1,:) = [a,b]; %#ok<AGROW>
        end
    end
end
