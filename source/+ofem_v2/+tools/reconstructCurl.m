function [X,U] = reconstructCurl(phys,u)
%RECONSTRUCT Reconstructs the solution from edge-elements
%   Output:
%   X = Coordinates of the evaluated points
%   U = Solution at the evaluation points
%   Input:
%   phys = ofem_v2.Physical_Problem type FE-Simulation data
%   u = FE-solution of the system
%   p = Order of approximation

    % First compute the global evaluation points
    [~,l] = ofem_v2.tools.gaussSimplex(phys.element.dim,1);%phys.element.degreeMass);
    Dk = phys.geometry.Dk;
    detD = phys.geometry.detD;
    Ne = phys.geometry.Nint;

        if phys.element.dim == 2
            l(3,:) = 1-sum(l,1);
        elseif phys.element.dim ==3    
            l(4,:) = 1-sum(l,1);
        end

    Nl = size(l,1);
    Np = size(l,2);
    elco = reshape(phys.geometry.co(:,:,phys.geometry.el(:,:)'),[],Nl,Ne);
    X = zeros(3,Ne,Np);
    for i = 1:Np
        X(:,:,i) = reshape(elco*l(:,i),3,[]);
    end    

    elu = ofem_v2.tools.matrixarray(reshape(u(phys.DOFs.el2DOF(:,:)'),[],1,Ne));
    for q = 1:Np
        dphi(:,:,1) = phys.element.curlN{1}(l(1,q),l(2,q),l(3,q));
        dphi(:,:,2) = phys.element.curlN{2}(l(1,q),l(2,q),l(3,q));
        dphi = ofem_v2.tools.matrixarray(dphi(:,:,phys.geometry.refTet));
        U(:,:,q) = reshape((Dk*dphi)*elu*ofem_v2.tools.matrixarray(1/(detD)),3,[]);
    end
    
    
end

