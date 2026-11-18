function [mur,R,d] = addShells(mesh,varargin)
%ADDSHELLS materials of the improvised asymptotic boundary condition
%(Meeker 2013) for the A formulation (H(curl), stiffness 1/mu)
%
%   [mur,R,d] = ofem_v2.boundary.iabc.addShells(mesh)
%   [mur,R,d] = ofem_v2.boundary.iabc.addShells(mesh, 'Layers', {'Layer1','Layer2','Layer3'}, 'Mass', 1)
%   [mur,R,d] = ofem_v2.boundary.iabc.addShells(mesh, 'Mur', [1.01053 0.709899 6.91301])
%
%   The shells (parts 'Layer1' ... 'LayerN' by default, N = order of the
%   ABC) are found in the mesh, R and d are measured from it and the
%   permeabilities mur are computed for delta = d/R. Each shell gets a
%   material with stiff = 1/(mur*mu0) and mass = 'Mass' (regularisation,
%   as for the other materials of the problem). The outer boundary of the
%   last shell needs n x A = 0 (DirichletEdge). 'Mur' prescribes the
%   permeabilities instead (e.g. to reproduce older results).
    ip = inputParser;
    ip.addParameter('Layers',[]);
    ip.addParameter('Mass',1);
    ip.addParameter('Mur',[]);
    ip.parse(varargin{:});
    layers = ip.Results.Layers;
    if isempty(layers)
        layers = ofem_v2.boundary.iabc.layerNames(mesh);
    end
    [R,d] = ofem_v2.boundary.iabc.shellGeometry(mesh,layers);
    if isempty(ip.Results.Mur)
        mur = ofem_v2.boundary.iabc.shellPermeability(d/R,numel(layers));
    else
        mur = ip.Results.Mur(:);
    end
    for k = 1:numel(layers)
        m = ofem_v2.materials.Material;
        m.stiff = 1/(mur(k)*m.mu0);
        m.mass = ip.Results.Mass;
        mesh.setMaterial(layers{k},m);
    end
end
