%% Eddy currents in a copper wire: frequency sweep with iterative solvers
%
% Time-harmonic eddy current problem for the magnetic vector potential A
%
%   curl(nu curl A) + 1i*omega*kappa*A = J0     in the box
%                                n x A = 0      on the outer boundary
%
% with an impressed current density J0 = (0,0,J0) in the wire. The air gets
% a small conductivity (kappa = 1) as regularisation, as is common for
% A-formulations. Discretised with lowest-order H(curl) elements:
%
%   (S + 1i*omega*M) u = b,   S: nu-weighted curl-curl, M: kappa-weighted mass
%
% The system is solved with BiCGStab and two preconditioners:
%
% AdditiveJacobi (frequency sweep mode)
%   element-block Jacobi of S + omega*M plus a correction on the gradients.
%   The blocks of S and M are extracted once, setFrequency only inverts
%   them. Since G'*S*G = 0, the gradient matrix G'*A*G = 1i*omega*G'*M*G is
%   factorised once for all frequencies (incomplete Cholesky).
% HiptmairXu
%   auxiliary space preconditioner with algebraic multigrid on the nodal
%   (H1) auxiliary problems; its iteration count does not grow with the
%   mesh size. Built for S + omega*M at every frequency.
%
% Before solving, the gradient part of the right hand side is removed
% (AdditiveJacobi.gradientSolve), which keeps the discrete current source
% divergence free.
%
% Outputs per frequency: iterations and times of both solvers, the
% difference of the solutions, the time averaged magnetic energy
% Wm = 1/4 u'*S*u and the eddy current losses in the wire
% P = 1/2 omega^2 u'*M_wire*u. With an impressed current density the
% induced current cancels J0 at high frequencies (the field is pushed
% out), so Wm decreases and P approaches 1/2*J0^2*V_wire/kappa.

close all
clear

h     = 1;                          % mesh refinement
f     = [1e2 1e3 1e4 1e5 1e6];      % frequencies [Hz]
J0    = 1e6;                        % impressed current density [A/m^2]
tol   = 1e-8;                       % BiCGStab tolerance

%% Mesh
file = 'wire';
system(sprintf('gmsh %s.geo -setnumber h %g -3 -v 0 -o %s.msh',file,h,file));
mesh = ofem_v2.Geometry;
mesh.load_from_msh(file);
mesh.reorderAC;
mesh.create_faces;
mesh.create_edges;
mesh.connectFa2Ed;
mesh.jacobiandata;

%% Materials: mu holds the reluctivity 1/mu, kappa the conductivity
copper = ofem_v2.materials.Material;
copper.mu = 1/copper.mu0;
copper.kappa = 5.8e7;
air = ofem_v2.materials.Material;
air.mu = 1/air.mu0;
air.kappa = 1;                      % regularisation
mesh.setMaterial('Wire',copper);
mesh.setMaterial('Air',air);

%% Problem
fe = ofem_v2.elements.loadFE('HCurl_3D_Order_0');
dofs = ofem_v2.DOFHandler(mesh);
dofs.attach(fe);

phys = ofem_v2.Physical_Problem(fe,mesh,1,1,0);
phys.setParaS('mu');
phys.setParaM('kappa');
phys.setBoundaryCondition(ofem_v2.boundary.DirichletEdge([0;0;0],'PEC',mesh));
ofem_v2.Volume_force([0;0;J0],'Wire',mesh);
phys.attachDOFHandler(dofs);
dofs.generateDOFs;
phys.assemble(1);

S = phys.S;
M = phys.M;
b = phys.b;
free = dofs.freeDOFs;
Mwire = phys.MPart{strcmp(mesh.parts(1,:),'Wire')};
fprintf('%d tetrahedra, %d free DOFs\n',size(mesh.el,1),numel(free));

%% AdditiveJacobi: frequency independent setup
t = tic;
G = fe.AMS(phys);
aj = ofem_v2.preconditioners.AdditiveJacobi(dofs,S,M);
aj.gradientCorrection(G,M,free,dofs.AMSDOFs(phys),'Solver','ict');
fprintf('AdditiveJacobi setup (once): %.1fs\n\n',toc(t));

fprintf('%9s | %-24s | %-24s | %9s | %11s %11s\n','f [Hz]', ...
    'AdditiveJacobi: it, time','HiptmairXu: setup+solve','diff','Wm [J]','P [W]');
Wm = zeros(size(f));
P = zeros(size(f));
for k = 1:numel(f)
    omega = 2*pi*f(k);
    A = S + 1i*omega*M;
    Af = A(free,free);

    % remove the gradient part of the right hand side
    aj.setFrequency(omega);
    phi = aj.gradientSolve(aj.G'*b(free));
    bc = b(free) - 1i*omega*M(free,free)*(aj.G*phi);

    t = tic;
    [u1,fl1,~,it1] = bicgstab(Af,bc,tol,3000,@(r) aj.preconditioner(r));
    t1 = toc(t);

    t = tic;
    hx = ofem_v2.preconditioners.HiptmairXu.fromProblem(phys,S + omega*M);
    ts = toc(t);
    t = tic;
    [u2,fl2,~,it2] = bicgstab(Af,bc,tol,3000,@(r) hx.apply(r));
    t2 = toc(t);

    u = zeros(dofs.Nd,1);
    u(free) = u2;
    Wm(k) = real(u'*S*u)/4;
    P(k) = omega^2*real(u'*Mwire*u)/2;
    fprintf('%9.0e | %6.1f it %6.1fs%s       | %4.1f it %4.1f+%4.1fs%s      | %9.1e | %11.4e %11.4e\n', ...
        f(k),it1,t1,flagText(fl1),it2,ts,t2,flagText(fl2),norm(u1-u2)/norm(u2),Wm(k),P(k));
end

% check: at high frequencies the induced current cancels J0 in the wire,
% so the losses approach P = 1/2*J0^2*V_wire/kappa
Vwire = sum(abs(mesh.detD(mesh.parts{3,strcmp(mesh.parts(1,:),'Wire')})))/6;
fprintf('\nhigh frequency limit of the losses: %.4e W (computed at %.0e Hz: %.4e W)\n', ...
    J0^2*Vwire/(2*copper.kappa),f(end),P(end));

%% Plot
figure('Name','Eddy currents in a wire');
subplot(1,2,1); loglog(f,Wm,'o-'); grid on
xlabel('f [Hz]'); ylabel('W_m [J]'); title('magnetic energy')
subplot(1,2,2); loglog(f,P,'o-'); grid on
xlabel('f [Hz]'); ylabel('P [W]'); title('eddy current losses in the wire')

function s = flagText(flag)
    if flag == 0, s = ' '; else, s = '!'; end
end
