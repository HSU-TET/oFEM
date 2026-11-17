function [err,nd] = rect(mesh,order)

	
	%fe = ofem_v2.elements.loadFE('HCurl_3D_Order_3');
	fe = ofem_v2.elements.HCurlElement(2,order);
	fe.computeBasis;

	
	dofs = ofem_v2.DOFHandler(mesh);
	dofs.attach(fe);
	
	PEC = ofem_v2.boundary.DirichletEdge([0;0],'PEC',mesh);
	
	resonator = ofem_v2.Physical_Problem(fe,mesh,1,1,0);
	
	resonator.setParaS('stiff');
	resonator.setParaM('mass');
	
	resonator.setBoundaryCondition(PEC);
	
	resonator.attachDOFHandler(dofs);
	dofs.generateDOFs;
	
	resonator.assemble;
	
	Y_g = fe.AMS(resonator);
	amsDOFS = dofs.AMSDOFs(resonator);
	
	S = resonator.S;
	M = resonator.M;
	
	x = resonator.u;
	
	X = randn(length(x),14);
	x = zeros(size(X));
	
	Y_g = Y_g(:,amsDOFS);
	S_g = Y_g'*M*Y_g;
	
	[x(dofs.freeDOFs,:),l,~,~] = ofem_v2.solvers.lobpcgModProj(X,S,M,[],7,50,1e-9,Y_g,S_g,dofs.freeDOFs);
	exact = [1,1,2,4,4,5,5,8,9,9,10,10,13,13]';
	err = l-exact;
	nd = dofs.Nd;
end

