function B = traceBasis(td, geom, fun, map)
%TRACEBASIS volume basis of the adjacent elements at the quadrature points
%of traceData
%
%   B = ofem_v2.boundary.traceBasis(td, geom, fun, map)
%
%   fun  basis on the two reference elements, e.g. fe.phi (H1) or fe.N
%        (HCurl / HDiv): fun{v}(u,v[,w]) returns nc x Ns
%   map  'none'          H1 values (default)
%        'covariant'     DinvT*N   (H(curl))
%        'contravariant' Dk*N/detD (H(div))
%
%   B  nc x Ns x Nq x Nf, the DOFs are el2DOF(td.elem,:)

    if nargin < 4
        map = 'none';
    end
    [~,Nq,nl] = size(td.P);
    c = num2cell(td.P(:,1,1));
    [nc,Ns] = size(fun{1}(c{:}));
    % evaluate once per local face and variant
    T = zeros(nc,Ns,Nq,nl,2);
    for v = 1:2
        for k = 1:nl
            for q = 1:Nq
                c = num2cell(td.P(:,q,k));
                T(:,:,q,k,v) = fun{v}(c{:});
            end
        end
    end
    e = td.elem(:);
    T = reshape(T,nc,Ns,Nq,nl*2);
    B = T(:,:,:,td.face(:) + nl*(double(geom.refTet(e))-1));     % nc x Ns x Nq x Nf
    switch map
        case 'none'
        case 'covariant'
            B = mapEach(double(geom.DinvT(:,:,e)),B);
        case 'contravariant'
            B = mapEach(double(geom.Dk(:,:,e))./reshape(double(geom.detD(:,:,e)),1,1,[]),B);
        otherwise
            error('ofem_v2:boundary:traceBasis','Unknown map ''%s''.',map);
    end
end

function B = mapEach(A,B)
    % B(:,:,q,f) = A(:,:,f)*B(:,:,q,f)
    [~,Ns,Nq,Nf] = size(B);
    B = permute(B,[1,2,4,3]);                       % nc x Ns x Nf x Nq
    out = zeros(size(A,1),Ns,Nf,Nq);
    for q = 1:Nq
        out(:,:,:,q) = pagemtimes(A,B(:,:,:,q));
    end
    B = permute(out,[1,2,4,3]);
end
