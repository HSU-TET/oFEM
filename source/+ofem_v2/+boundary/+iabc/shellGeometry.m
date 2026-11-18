function [R,d,c] = shellGeometry(mesh,layers)
%SHELLGEOMETRY radius R of the inner domain, shell thickness d and centre
%c of the IABC shells, read from the mesh
%
%   [R,d,c] = ofem_v2.boundary.iabc.shellGeometry(mesh, {'Layer1','Layer2','Layer3'})
%
%   The shells have to be concentric spherical shells of equal thickness
%   (Meeker's construction); a warning is issued otherwise.
    co = reshape(mesh.co,3,[]);
    nodesOf = @(name) unique(mesh.el(mesh.parts{3,partIndex(mesh,name)},1:4));
    outer = co(:,nodesOf(layers{end}));
    c = (min(outer,[],2) + max(outer,[],2))/2;          % centre of the outer sphere
    N = numel(layers);
    rIn = zeros(1,N); rOut = zeros(1,N);
    for k = 1:N
        r = vecnorm(co(:,nodesOf(layers{k})) - c);
        rIn(k) = min(r); rOut(k) = max(r);
    end
    R = rIn(1);
    d = (rOut(end) - R)/N;
    expected = R + (0:N)*d;
    if max(abs([rIn, rOut(end)] - expected([1:N, N+1]))) > 0.05*d
        warning('iabc:shells', ['The IABC shells are not of equal thickness ' ...
            '(radii %s), the permeabilities assume d = %g.'],mat2str([rIn rOut(end)],6),d);
    end
end

function i = partIndex(mesh,name)
    i = find(strcmp(mesh.parts(1,:),name),1);
    if isempty(i), error('iabc:part','The mesh has no part "%s".',name); end
end
