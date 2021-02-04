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
    l(4,:) = 1-sum(l,1);
    Nl = size(l,1);
    Np = size(l,2);
    elco = reshape(phys.geometry.co(:,:,phys.geometry.el(:,:)'),[],Nl,Ne);
    X = zeros(3,Ne,Np);
    for i = 1:Np
        X(:,:,i) = reshape(elco*l(:,i),3,[]);
    end
    elu = ofem_v2.tools.matrixarray(reshape(u(phys.DOFs.el2DOF(:,:)'),[],1,Ne));
    for q = 1:Np
        phi(:,:,1) = phys.element.N{1}(l(1,q),l(2,q),l(3,q));
        phi(:,:,2) = phys.element.N{2}(l(1,q),l(2,q),l(3,q));
        phi = ofem_v2.tools.matrixarray(phi(:,:,phys.geometry.refTet));
        U(:,:,q) = reshape((DinvT*phi)*elu,3,[]);
    end
        
end

