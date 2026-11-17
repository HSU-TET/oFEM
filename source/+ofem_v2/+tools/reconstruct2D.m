function [X,U] = reconstruct2D(phys,u,order)
%RECONSTRUCT2D Reconstructs the solution from 2D edge elements
%   Same as ofem_v2.tools.reconstruct, which handles 2D and 3D (this
%   function used the matrix product * on 3D arrays and could not run).
%   Output:
%   X = Coordinates of the evaluated points (2 x Ne x Np)
%   U = Solution at the evaluation points (2 x Ne x Np)
    [X,U] = ofem_v2.tools.reconstruct(phys,u,order);
end
