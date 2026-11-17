function G = geometryFactor(T,mat,fac)
%GEOMETRYFACTOR per element factor G = T'*mat*T.*fac for
%refTensorMatrices.
%
%   T    transformation of the reference basis (d x nc x Ne), [] for none
%   mat  material: scalar, matrix (nc x nc) or per element (nc x nc x Ne)
%   fac  per element factor (1 x 1 x Ne), e.g. |detD| or 1/|detD|
    if isempty(T)
        G = mat.*fac;
        return
    end
    if isscalar(mat)
        G = mat*pagemtimes(T,'transpose',T,'none');
    else
        G = pagemtimes(pagemtimes(T,'transpose',mat,'none'),T);
    end
    G = G.*fac;
end
