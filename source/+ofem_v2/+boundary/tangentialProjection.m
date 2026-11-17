function u = tangentialProjection(phys, cells, edges, faces, value)
%TANGENTIALPROJECTION DOF values of a prescribed tangential trace (H(curl))
%
%   u = ofem_v2.boundary.tangentialProjection(phys, cells, edges, faces, value)
%
%   Solves the L2 projection of the tangential trace of value onto the
%   trace of the DOFs of the boundary edges (and faces in 3D),
%
%       sum_j (n x phi_j, n x phi_i)_Gamma c_j = (n x g, n x phi_i)_Gamma,
%
%   which is exact for data in the discrete trace space (e.g. constant
%   fields) and works for any degree in 2D and 3D. The other DOFs of u
%   are zero.
%
%   value  constant vector (dim x 1), or a function handle g(X) or
%          g(X,n) of the points X (dim x 1 x N) and the outer unit
%          normals n (dim x 1 x N), returning dim x 1 x N

    D = phys.DOFs;
    g = phys.geometry;
    fe = phys.element;
    u = sparse(D.Nd,1);
    bd = [];
    if ~isempty(D.e2DOF) && ~isempty(edges)
        bd = [bd; reshape(D.e2DOF(edges,:),[],1)];
    end
    if ~isempty(D.f2DOF) && ~isempty(faces)
        bd = [bd; reshape(D.f2DOF(faces,:),[],1)];
    end
    bd = unique(bd);
    isFun = isa(value,'function_handle');
    if isempty(bd) || (~isFun && ~any(value(:)))
        return
    end

    order = 2*max(fe.degree,1) + 2*isFun;
    td = ofem_v2.boundary.traceData(g, cells, order);
    B = ofem_v2.boundary.traceBasis(td, g, fe.N, 'covariant');     % dim x Ns x Nq x Nf
    [dim,Ns,Nq,Nf] = size(B);

    % data at the quadrature points: dim x 1 x Nq x Nf
    G = zeros(dim,1,Nq,Nf);
    for q = 1:Nq
        if isFun
            if nargin(value) >= 2
                v = value(td.X(:,q,:),td.normal);
            else
                v = value(td.X(:,q,:));
            end
            G(:,1,q,:) = reshape(double(v),dim,1,1,Nf);
        else
            G(:,1,q,:) = repmat(value(:),1,1,1,Nf);
        end
    end

    % tangential parts: 2D phi.t with t = (-n2, n1), 3D phi - (phi.n) n
    n = reshape(td.normal,dim,1,1,Nf);
    if dim == 2
        t = [-n(2,:,:,:); n(1,:,:,:)];
        T = sum(t.*B,1);                    % 1 x Ns x Nq x Nf
        G = sum(t.*G,1);
    else
        T = B - n.*sum(n.*B,1);
        G = G - n.*sum(n.*G,1);
    end
    nt = size(T,1);

    % local matrices with the (positive) weights w(q)*jac
    sw = sqrt(reshape(td.w,1,1,Nq).*reshape(td.jac,1,1,1,Nf));
    Tm = reshape(permute(T.*sw,[1,3,2,4]),nt*Nq,Ns,Nf);
    Gm = reshape(permute(G.*sw,[1,3,2,4]),nt*Nq,1,Nf);
    Me = pagemtimes(Tm,'transpose',Tm,'none');
    Fe = pagemtimes(Tm,'transpose',Gm,'none');

    dofs = D.el2DOF(td.elem,:);
    I = repmat(dofs,1,Ns)';
    J = repelem(dofs',Ns,1);
    M = sparse(I(:),J(:),Me(:),D.Nd,D.Nd);
    b = sparse(reshape(dofs',[],1),1,Fe(:),D.Nd,1);
    u(bd) = M(bd,bd)\b(bd);
end
