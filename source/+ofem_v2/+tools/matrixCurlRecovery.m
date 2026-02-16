function [Mx,My,Mz,M] = matrixCurlRecovery(mesh,feCurl,feH)
    refTet = mesh.refTet;
    [w,l] = ofem_v2.tools.gaussSimplex(feCurl.dim,feCurl.degreeMass);
    dofsC = mesh.el2ed;
    dofsN = mesh.el;
    detD = mesh.detD;
    Dk = mesh.Dk;

    nN = size(dofsN,2);
    nC = size(dofsC,2);
    Ne = mesh.Nint;
    Nq = length(w);
    
    Mx = zeros(nN,nC,Ne);
    My = zeros(nN,nC,Ne);
    Mz = zeros(nN,nC,Ne);

    M = zeros(nN,nN,Ne);
    
    for q=1:Nq
        cnt = ones(size(l(:,q),1),1);
        lTemp = mat2cell(l(:,q),cnt);
        phiC(:,:,1) = feCurl.curlN{1}(lTemp{:});
        phiC(:,:,2) = feCurl.curlN{2}(lTemp{:});
        phiC = pagemtimes(Dk,phiC(:,:,refTet));
        cnt = ones(size(l(:,q),1),1);
        lTemp = mat2cell(l(:,q),cnt);
        phiN(:,:,1) = feH.phi{1}(lTemp{:});
        phiN(:,:,2) = feH.phi{2}(lTemp{:});
        phiN = phiN(:,:,refTet);
        Mx = Mx + w(q)*pagemtimes(phiN,'transpose',phiC(1,:,:),'none');
        My = My + w(q)*pagemtimes(phiN,'transpose',phiC(2,:,:),'none');
        Mz = Mz + w(q)*pagemtimes(phiN,'transpose',phiC(3,:,:),'none');
        M = M + w(q)*pagemtimes(phiN,'transpose',phiN,'none');
    end

    % Mx = pagemtimes(Mx,abs(detD));
    % My = pagemtimes(My,abs(detD));
    % Mz = pagemtimes(Mz,abs(detD));
    M = pagemtimes(M,abs(detD));
    
    I = repmat(dofsN,1,size(Mx,2))';
    %I = I(:);
    J = repelem(dofsC',size(Mx,1),1);
    IM = repmat(dofsN,1,size(M,1))';
    JM = repelem(dofsN',size(M,1),1);
    %J = J(:);
    
    Mx = sparse(I(:),J(:),Mx(:),mesh.Nco,mesh.Ned);
    My = sparse(I(:),J(:),My(:),mesh.Nco,mesh.Ned);
    Mz = sparse(I(:),J(:),Mz(:),mesh.Nco,mesh.Ned);

    M = sparse(IM(:),JM(:),M(:),mesh.Nco,mesh.Nco);
    
    % 			figure
    %             spy(M)
end

