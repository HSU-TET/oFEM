function [X,U] = reconstruct(phys,u,order)
%RECONSTRUCT Reconstructs the solution from edge-elements
%   Output:
%   X = Coordinates of the evaluated points
%   U = Solution at the evaluation points
%   Input:
%   phys = ofem_v2.Physical_Problem type FE-Simulation data
%   u = FE-solution of the system
%   order = Order of approximation

    % First compute the global evaluation points
    [~,l] = ofem_v2.tools.gaussSimplex(phys.element.dim,order);%phys.element.degreeMass);
    DinvT = phys.geometry.DinvT;
    Ne = phys.geometry.Nint;
    dim = phys.geometry.dim;
    if dim == 3
        l(4,:) = 1-sum(l,1);
    end
    if dim == 2
        l(3,:) = 1-sum(l,1);
    end
    Nl = size(l,1);
    Np = size(l,2);
    elco = reshape(phys.geometry.co(:,:,phys.geometry.el(:,:)'),[],Nl,Ne);
    X = zeros(dim,Ne,Np);
    for i = 1:Np
        X(:,:,i) = reshape(pagemtimes(elco,l(:,i)),dim,[]);
    end
    elu = reshape(u(phys.DOFs.el2DOF(:,:)'),[],1,Ne);
    for q = 1:Np
        cnt = ones(size(l(1:end-1,q),1),1);
        lTemp = mat2cell(l(1:end-1,q),cnt);
        phi(:,:,1) = phys.element.N{1}(lTemp{:});
        phi(:,:,2) = phys.element.N{2}(lTemp{:});
        phi = pagemtimes(DinvT,phi(:,:,phys.geometry.refTet));
        U(:,:,q) = reshape(pagemtimes(pagemtimes(DinvT,phi),elu),dim,[]);
    end
        
end

