function [I,idx,areas] = integrateOverROI(mesh,ROI,phi,dir,fe,dofs)
	for i = 1:size(mesh.roi,2)
		if startsWith(mesh.roi{1,i},ROI)
			faces = mesh.roi{2,i};
		end
	end
	dir = ofem_v2.tools.matrixarray(repmat(dir,1,1,size(faces,1)));
	faces = sort(faces,2);
	v1 = mesh.co(:,:,faces(:,2))-mesh.co(:,:,faces(:,1));
	v2 = mesh.co(:,:,faces(:,3))-mesh.co(:,:,faces(:,1));
	areas = sqrt(dot(cross(v1,v2),cross(v1,v2)))/2;
	[~,fidx] = ismember(faces,mesh.fa,'rows');

	tmp = zeros(size(fidx,1),2);
	for i=1:length(fidx)
		tmp(i,:) = mod(find(mesh.el2fa==fidx(i)),mesh.Nint);
	end
	centers1 = (reshape(mesh.co(:,:,mesh.el(tmp(:,1),:)'),[],4,length(fidx))*[1/4;1/4;1/4;1/4])'*dir;
	centers2 = (reshape(mesh.co(:,:,mesh.el(tmp(:,2),:)'),[],4,length(fidx))*[1/4;1/4;1/4;1/4])'*dir;
	ref = (centers1(1,:,:)>centers2(1,:,:))+1;
	for i = 1:length(ref)
		idx(i) = tmp(i,ref(i));
	end
	kappa = mesh.parts{2,1}.kappa;
	[~,facenr] = ismember(fidx,mesh.el2fa(:));
	facenr = facenr/mesh.Nint;
	facenr = floor(facenr)+1;
	functs = [1,2,3;1,2,4;1,3,4;2,3,4];
	functs = functs(facenr,:);

	DinvT = mesh.DinvT(:,:,idx);
	uElem = phi(dofs.el2DOF(idx,:));
	uElem = ofem_v2.tools.matrixarray(reshape(uElem',size(uElem,2),1,[]));
	dphi(:,:,1) = fe.dPhi{1}(1/3,1/3,1/3);
	dphi = ofem_v2.tools.matrixarray(dphi);
	du = (DinvT*dphi)*uElem;

	E = -du;
	J = kappa*E;

	I = sum((J'*dir)*areas);
end



















