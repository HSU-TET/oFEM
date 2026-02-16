function u_rec = nodalCurlRecovery(mesh,u)
    % Recovers the curl of solution vector for edge elements at the nodes
    % Returns Nco X 3 Array containing the vectorial values 
    % Only works for constant mu
    feCurl = ofem_v2.elements.loadFE('HCurl_3D_Order_0');
    feH = ofem_v2.elements.loadFE('H1_3D_Order_1');
    
    [Mx,My,Mz,M] = ofem_v2.tools.matrixCurlRecovery(mesh,feCurl,feH);
    
    u_x = M\(Mx*u);
    u_y = M\(My*u);
    u_z = M\(Mz*u);
    
    u_rec = [u_x,u_y,u_z];
end