function td = traceData(geom, cells, order)
%TRACEDATA quadrature on boundary cells through their adjacent elements
%
%   td = ofem_v2.boundary.traceData(geom, cells, order)
%
%   cells  boundary cells as node lists (Nf x 2 edges in 2D, Nf x 3
%          triangles in 3D), e.g. geom.bd{2,i}
%   order  polynomial degree the face quadrature integrates exactly
%
%   For every cell the adjacent element and the local face (edge in 2D)
%   of the reference element are determined. The quadrature rule of the
%   reference face is mapped into the element, so the volume basis is
%   evaluated at the same reference points for all cells with the same
%   local face and refTet variant (see traceBasis). This gives the trace
%   of the volume basis with the global DOFs el2DOF(td.elem,:) and the
%   right orientation for any polynomial degree.
%
%   td.w       weights of the reference face rule (Nq x 1)
%   td.X       physical quadrature points (dim x Nq x Nf)
%   td.elem    adjacent element of each cell (Nf x 1)
%   td.face    local face / edge of the cell in that element (Nf x 1)
%   td.P       reference points per local face (dim x Nq x nfaces)
%   td.normal  unit normal pointing out of the adjacent element (dim x 1 x Nf)
%   td.jac     |face| / |reference face| (1 x 1 x Nf), so that
%              int_cell f = sum_q w(q) f(X(:,q,f)) jac(f)
%   td.meas    length / area of the cells (1 x 1 x Nf)

    dim = geom.dim;
    cells = sort(cells,2);
    switch dim
        case 2
            [~,ent] = ismember(cells,geom.ed,'rows');
            el2x = geom.el2ed;
            loc = [1,2;1,3;2,3];        % local edges of create_edges
        case 3
            [~,ent] = ismember(cells,geom.fa,'rows');
            el2x = geom.el2fa;
            loc = [1,2,3;1,2,4;1,3,4;2,3,4];   % local faces of create_faces
        otherwise
            error('ofem_v2:boundary:traceData','Only 2D and 3D meshes are supported.');
    end
    if any(ent == 0)
        error('ofem_v2:boundary:traceData','Boundary cells are not edges / faces of the mesh.');
    end

    % adjacent element and its local face of every cell
    Ne = size(el2x,1);
    nl = size(el2x,2);
    owner = zeros(max(el2x(:)),1);
    lface = owner;
    owner(el2x(:)) = repmat((1:Ne)',nl,1);
    lface(el2x(:)) = repelem((1:nl)',Ne,1);
    td.elem = owner(ent);
    td.face = lface(ent);

    % face rule mapped to every local face of the reference element
    [w,s] = ofem_v2.tools.gaussSimplex(dim-1,order);
    if dim == 2
        s = s(:)';                      % 1 x Nq
    end
    td.w = w(:);
    Nq = numel(td.w);
    Vref = [zeros(dim,1), eye(dim)];    % vertex 1 at the origin
    td.P = zeros(dim,Nq,nl);
    for k = 1:nl
        a = Vref(:,loc(k,1));
        td.P(:,:,k) = a + (Vref(:,loc(k,2:end)) - a)*s;
    end

    % physical points
    e = td.elem;
    Nf = numel(e);
    x1 = geom.co(:,1,geom.el(e,1));
    td.X = double(x1) + pagemtimes(double(geom.Dk(:,:,e)),td.P(:,:,td.face));

    % normals and face measures from the cell geometry
    c = double(geom.co(:,1,cells'));
    c = reshape(c,dim,dim,Nf);          % vertices of each cell
    t1 = c(:,2,:) - c(:,1,:);
    if dim == 2
        n = [t1(2,:,:); -t1(1,:,:)];
        nrm = sqrt(sum(n.^2,1));
        td.jac = nrm;
        td.meas = nrm;
    else
        n = cross(t1,c(:,3,:) - c(:,1,:),1);
        nrm = sqrt(sum(n.^2,1));
        td.jac = nrm;                   % reference triangle has area 1/2
        td.meas = nrm/2;
    end
    n = n./nrm;
    centroid = mean(reshape(double(geom.co(:,1,geom.el(e,:)')),dim,dim+1,Nf),2);
    sgn = sign(sum(n.*(c(:,1,:) - centroid),1));
    sgn(sgn == 0) = 1;
    td.normal = n.*sgn;
end
