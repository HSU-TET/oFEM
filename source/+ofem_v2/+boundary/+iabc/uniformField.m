function uniformField(mesh,B0,varargin)
%UNIFORMFIELD uniform background field B0 inside the IABC shells
%
%   ofem_v2.boundary.iabc.uniformField(mesh, B0)
%   ofem_v2.boundary.iabc.uniformField(mesh, B0, {'Layer1','Layer2','Layer3'})
%
%   A uniform coercivity Hc = -3/2 B0/mu0 in all shells acts as a surface
%   current Hc x n on the inner sphere r = R (the contributions at the
%   shell interfaces cancel, the outer one is removed by n x A = 0). With
%   the exterior emulated by the IABC (exact for the dipole harmonic) the
%   field inside the sphere is the field of a uniformly magnetised sphere,
%   2/3 mu0 |Hc| = |B0|.
    if nargin < 3
        layers = ofem_v2.boundary.iabc.layerNames(mesh);
    else
        layers = varargin{1};
    end
    mu0 = 4e-7*pi;
    for k = 1:numel(layers)
        mesh.setForce(layers{k},ofem_v2.Coercivity(-1.5*B0(:)/mu0,layers{k},mesh));
    end
end
