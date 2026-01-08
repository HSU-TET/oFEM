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
    dim = phys.geometry.dim;

    if dim == 2
        l(3,:) = 1-sum(l,1);
        outsize = 1;
    elseif dim ==3    
        l(4,:) = 1-sum(l,1);
        outsize = 3;
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
        dphi(:,:,1) = phys.element.curlN{1}(lTemp{:});
        dphi(:,:,2) = phys.element.curlN{2}(lTemp{:});
        dphi = dphi(:,:,phys.geometry.refTet);
        if size(dphi,1)==3
            dphi = pagemtimes(Dk,dphi);
        end
        U(:,:,q) = reshape(pagemtimes(pagemtimes(dphi,elu),pageinv(detD)),outsize,[]);
        
    end
    
    
end

