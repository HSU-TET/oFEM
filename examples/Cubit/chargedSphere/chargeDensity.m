function rhs = chargeDensity(X)
    rho = ofem_v2.tools.matrixarray(sqrt(dot(X,X,1)));

    rhs = 0*rho;
    rhs(:,:,rho(:)<=8)=1;
%     rhs(:,:,rho(:)<=8) = 100.*exp(-dot(X,X,1));
end