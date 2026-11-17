function Me = refTensorMatrices(basis,w,l,refTet,G)
%REFTENSORMATRICES element matrices on affine elements via reference
%tensors.
%
%   Me = refTensorMatrices(basis,w,l,refTet,G) returns the Ns x Ns x Ne
%   element matrices
%
%       Me(:,:,e) = sum_q w(q) * B_q' * G(:,:,e) * B_q,
%
%   where B_q = basis{refTet(e)}(l(:,q)) is the nc x Ns matrix of the
%   reference basis (or its derivatives) at the quadrature point l(:,q).
%   G (nc x nc x Ne) contains everything that is constant per element:
%   the transformation, the material and the Jacobian determinant, e.g.
%   G = DinvT'*mat*DinvT*|detD| for a covariant mass matrix.
%
%   Since G does not depend on the quadrature point, the matrices are
%   Me(:,:,e) = sum_ab G(a,b,e) * R_ab with the reference tensors
%   R_ab = sum_q w(q) * B_q(a,:)' * B_q(b,:), so the work per element is
%   a single matrix product instead of a loop over quadrature points.

    nc = size(G,1);
    Ne = size(G,3);
    nv = numel(basis);
    c = num2cell(l(:,1));
    Ns = size(basis{1}(c{:}),2);
    R = zeros(Ns*Ns,nc*nc,nv);
    for v = 1:nv
        for q = 1:numel(w)
            c = num2cell(l(:,q));
            B = basis{v}(c{:});
            for b = 1:nc
                for a = 1:nc
                    A = B(a,:)'*B(b,:);
                    R(:,a+(b-1)*nc,v) = R(:,a+(b-1)*nc,v) + w(q)*A(:);
                end
            end
        end
    end
    G = reshape(G,nc*nc,Ne);
    Me = zeros(Ns*Ns,Ne);
    for v = 1:nv
        e = refTet == v;
        if any(e)
            Me(:,e) = R(:,:,v)*G(:,e);
        end
    end
    Me = reshape(Me,Ns,Ns,Ne);
end
