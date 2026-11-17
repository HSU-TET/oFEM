%% Pulsating flow past an off-center cylinder with tracer transport
%
% 1) Velocity: potential (Darcy) flow as in ../../Static/CylinderFlow,
%
%       u + grad p = f(t),  div u = 0,  u.n = 0 on walls and cylinder,
%       p = 0 at inlet and outlet,
%
%    driven by the pulsating force f(t) = G(t)*(1,0) with
%    G(t) = G0*(1 + 0.5*sin(2*pi*t/Tp)). Darcy flow has no inertia, so
%    u(t) = G(t)/G0 * u0 with u0 the solution for G0: the H(div) problem
%    is solved once. The cylinder sits below the channel center line.
%
% 2) Tracer: a concentration c(x,t) transported by the flow,
%
%       dc/dt + div(u c) = 0,   c = c_in(y,t) at the inlet,
%
%    discretised with a piecewise constant upwind finite volume scheme on
%    the triangles and SSP-RK2 in time. The face fluxes are the exact
%    edge fluxes int_e u_h.n of the H(div) velocity. Because div u_h = 0
%    on every element, the fluxes of every triangle sum to zero: the
%    scheme conserves mass exactly and keeps a constant concentration
%    constant. Dye stripes are injected at the inlet for t < Tin.

close all
clear

order = 2;       % H(div) order of the velocity
h     = 2;       % mesh refinement
yc    = 0.13;    % cylinder center height (channel center: 0.205)
G0    = 1;       % mean driving force
Tp    = 0.5;     % period of the pulsation
Tin   = 0.6;     % dye is injected for t < Tin
T     = 3;       % final time
CFL   = 0.4;

% UCD export (ParaView) of concentration and velocity per triangle.
% exportEvery = 1 writes every time step: with h = 2 about 17000 files of
% 1.1 MB (~19 GB, ~10 min). Increase it to write every n-th step only.
exportEvery = 50;

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
file = fullfile(here,'channelOffset');
exportDir = fullfile(here,'export');

%% Mesh
[status,msg] = system(sprintf('gmsh "%s" -setnumber h %g -setnumber yc %g -2 -v 0 -o "%s.msh"', ...
    fullfile(here,'channel.geo'),h,yc,file));
if status ~= 0
    error('gmsh failed: %s',msg);
end
mesh = ofem_v2.Geometry;
mesh.load_from_msh(file);
mesh.reorderAC;
mesh.create_edges;
mesh.jacobiandata;

fluid = ofem_v2.materials.Material;
fluid.mass = 1;
mesh.setMaterial('Domain',fluid);

%% Steady velocity for G0
fe = ofem_v2.elements.loadFE(sprintf('HDiv_2D_Order_%d',order));
dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(fe);
phys = ofem_v2.Physical_Problem(fe,mesh,0,1,0);
phys.setParaM('mass');
phys.setBoundaryCondition(ofem_v2.boundary.NormalTrace('Wall',mesh));
phys.setBoundaryCondition(ofem_v2.boundary.NormalTrace('Cylinder',mesh));
phys.attachDOFHandler(dofs);
dofs.generateDOFs;
phys.assemble;

M = phys.M;
F = fe.volumeForce(phys,[G0;0],1:mesh.Nint);
B = divergenceMatrix(fe,phys,max(order-1,0));
free = dofs.freeDOFs;
Nf = numel(free);
x = [M(free,free), -B(:,free)'; -B(:,free), sparse(size(B,1),size(B,1))] ...
    \ [F(free); zeros(size(B,1),1)];
u0 = zeros(dofs.Nd,1);
u0(free) = x(1:Nf);

% element mean velocity for G0 (for the export), u(t) = G(t)/G0 * u0
[wq,lq] = ofem_v2.tools.gaussSimplex(2,2*order);
[~,Uq] = fe.evaluate(phys,u0,lq);
Ue0 = 2*sum(Uq.*reshape(wq,1,1,[]),3)';             % Ne x 2

%% Edge fluxes of the H(div) velocity
% Fl(k,K): outward flux of u0 through local edge k of triangle K
Fl = edgeFluxes(fe,phys,u0);
Ne = mesh.Nint;
area = abs(reshape(mesh.detD,1,[]))/2;
netFlux = abs(sum(Fl,1))./max(abs(Fl),[],1);

% neighbours over every local edge (0 on the boundary)
ed = mesh.el2ed;
owner = accumarray(ed(:),repmat((1:Ne)',3,1),[mesh.Ned,1],@(k) {k});
nb = zeros(3,Ne);
for K = 1:Ne
    for k = 1:3
        o = owner{ed(K,k)};
        o = o(o ~= K);
        if ~isempty(o), nb(k,K) = o; end
    end
end

% upwind operator dc/dt = A*c + b(t) for u0
out = Fl > 0;
[kk,KK] = find(out);
f = Fl(out);
NB = nb(out);
inner = NB > 0;
A = sparse(KK,KK,-f./area(KK)',Ne,Ne) ...
  + sparse(NB(inner),KK(inner),f(inner)./area(NB(inner))',Ne,Ne);

% inflow edges: flux and midpoint height
co = reshape(mesh.co,2,[]);
[kin,Kin] = find(Fl < 0 & nb == 0);
fin = Fl(sub2ind(size(Fl),kin,Kin));
lv = [1 2;1 3;2 3];
yin = zeros(size(Kin));
for i = 1:numel(Kin)
    v = mesh.el(Kin(i),lv(kin(i),:));
    yin(i) = mean(co(2,v));
end
bIn = @(cin) accumarray(Kin,-fin.*cin./area(Kin)',[Ne,1]);

% concentration at the inlet: 3 dye stripes while t < Tin
Hc = 0.41;
cIn = @(y,t) double(t < Tin & mod(floor(y/(Hc/7)),2) == 1);
s = @(t) 1 + 0.5*sin(2*pi*t/Tp);                    % G(t)/G0

fprintf('H(div) order %d, %d triangles\n',order,Ne);
fprintf('  max |sum of edge fluxes| / max edge flux per triangle : %.1e\n',max(netFlux));
fprintf('  constant concentration preserved, max |dc/dt|        : %.1e\n', ...
    max(abs(A*ones(Ne,1) + bIn(ones(size(Kin))))));

%% Time stepping (SSP-RK2)
dt = CFL*min(area'./sum(max(Fl,0),1)')/max(s(0:0.01:Tp));
nt = ceil(T/dt);
dt = T/nt;
c = zeros(Ne,1);
mass = zeros(nt+1,1);
balance = zeros(nt+1,1);                            % inflow - outflow
snap = [0.4 1.0 1.6 2.4];
C = zeros(Ne,numel(snap));
isnap = 1;
[~,Kout] = find(Fl > 0 & nb == 0);                  % outflow edges
fout = Fl(Fl > 0 & nb == 0);
L = @(c,t) s(t)*(A*c + bIn(cIn(yin,t)));
net = @(c,t) s(t)*(sum(-fin.*cIn(yin,t)) - sum(fout.*c(Kout)));

% export: files cylinderTransport_<step>.inp, ParaView loads them as a
% time series; the times of the steps are written to ..._times.csv
if exist(exportDir,'dir')
    delete(fullfile(exportDir,'cylinderTransport_*.inp'));
end
exportStep = @(n,c) exportUCD(mesh,exportDir,n,c,[s(n*dt)*Ue0, zeros(Ne,1)]);
exportStep(0,c);
exported = 0;
for n = 1:nt
    t = (n-1)*dt;
    c1 = c + dt*L(c,t);
    balance(n+1) = balance(n) + dt/2*(net(c,t) + net(c1,t+dt));
    c = (c + c1 + dt*L(c1,t+dt))/2;
    mass(n+1) = area*c;
    if mod(n,exportEvery) == 0
        exportStep(n,c);
        exported(end+1) = n; %#ok<SAGROW>
    end
    while isnap <= numel(snap) && n*dt >= snap(isnap)-dt/2
        C(:,isnap) = c;
        isnap = isnap+1;
    end
end
writematrix([exported(:), exported(:)*dt],fullfile(exportDir,'cylinderTransport_times.csv'));
fprintf('  time steps                                            : %d (dt = %.2e)\n',nt,dt);
fprintf('  exported files                                        : %d in %s\n',numel(exported),exportDir);
fprintf('  mass balance error max |mass - (inflow - outflow)|    : %.1e (max mass %.3f)\n', ...
    max(abs(mass-balance)),max(mass));
fprintf('  concentration range                                   : [%.3f, %.3f]\n',min(C(:)),max(C(:)));

%% Plot
tri = mesh.el(:,1:3);
figure('Name','Tracer in pulsating flow past an off-center cylinder','Position',[100 50 1000 800]);
for i = 1:numel(snap)
    subplot(numel(snap)+1,1,i)
    patch('Faces',tri,'Vertices',co','FaceVertexCData',C(:,i), ...
        'FaceColor','flat','EdgeColor','none');
    axis equal tight; clim([0 1]); colorbar
    title(sprintf('tracer concentration, t = %.1f',snap(i)))
end
subplot(numel(snap)+1,1,numel(snap)+1)
tt = (0:nt)*dt;
yyaxis left; plot(tt,mass); ylabel('tracer mass')
yyaxis right; plot(tt,s(tt)); ylabel('G(t)/G_0')
xlabel('t'); grid on

%% Local functions

function exportUCD(mesh,folder,n,c,u) %#ok<INUSD> used inside evalc
    % concentration and velocity per triangle at time step n (the console
    % output of export_UCD is suppressed)
    evalc(['mesh.export_UCD(folder,sprintf(''cylinderTransport_%05d'',n),' ...
        '{''concentration'',c,''1'',''Cell''},{''velocity'',u,''m/s'',''Cell''})']);
end

function Fl = edgeFluxes(fe,phys,u)
    % outward fluxes int_e u.n through the 3 local edges of every element.
    % The Piola map preserves fluxes up to the orientation: flux =
    % sign(detD) * reference flux.
    g = phys.geometry;
    [wq,tq] = ofem_v2.tools.gaussSimplex(1,fe.degree+1);
    % reference edges [1 2], [1 3], [2 3]: points and length*outward normal
    pts = {@(t) [t;0], @(t) [0;t], @(t) [1-t;t]};
    nrm = {[0;-1], [-1;0], [1;1]};
    Ns = fe.DOFsPerElement;
    Phi = zeros(Ns,3,2);
    for var = 1:2
        for k = 1:3
            for q = 1:numel(wq)
                xq = pts{k}(tq(q));
                Phi(:,k,var) = Phi(:,k,var) + wq(q)*(fe.N{var}(xq(1),xq(2))'*nrm{k});
            end
        end
    end
    ue = reshape(full(u(phys.DOFs.el2DOF')),Ns,1,[]);
    Fl = pagemtimes(Phi(:,:,g.refTet),'transpose',ue,'none');
    Fl = reshape(Fl.*sign(g.detD),3,[]);
end

function B = divergenceMatrix(fe,phys,k)
    % B(i,j) = int q_i div(v_j), q_i monomials of degree <= k in the
    % reference coordinates of each element (discontinuous P_k)
    g = phys.geometry;
    Ne = g.Nint;
    E = zeros(0,2);
    for a = 0:k
        for b = 0:k-a
            E(end+1,:) = [a,b]; %#ok<AGROW>
        end
    end
    nq = size(E,1);
    [w,l] = ofem_v2.tools.gaussSimplex(2,max(k+max(fe.degree-1,0),1));
    Ns = fe.DOFsPerElement;
    Bl = zeros(nq,Ns,Ne);
    for q = 1:numel(w)
        dphi = cat(3,fe.divN{1}(l(1,q),l(2,q)),fe.divN{2}(l(1,q),l(2,q)));
        dphi = dphi(:,:,g.refTet);
        Bl = Bl + w(q)*prod(l(:,q)'.^E,2).*dphi;
    end
    Bl = Bl.*sign(g.detD);
    rows = reshape(1:nq*Ne,nq,Ne);
    I = repmat(reshape(rows,nq,1,Ne),1,Ns,1);
    J = repmat(reshape(phys.DOFs.el2DOF',1,Ns,Ne),nq,1,1);
    B = sparse(I(:),J(:),Bl(:),nq*Ne,phys.DOFs.Nd);
end
