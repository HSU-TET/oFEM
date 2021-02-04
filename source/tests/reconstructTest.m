i = 2;
%u = (phys.S-ksq(i)*phys.M)*v1(:,i);
%u = v1(:,i);
%[x,v] = ofem_v2.tools.reconstruct(phys,full(v1(:,i)));
[x,v] = ofem_v2.tools.reconstructCurl(phys,full(u));
thres = 1296;
% plot3(x(1,:),x(2,:),x(3,:),'kx')
% hold on
% tetramesh(TR,'FaceAlpha',0.2)
X = x(1,:,:);
X = X(:);
Y = x(2,:,:);
Y = Y(:);
Z = x(3,:,:);
Z = Z(:);
Ux = v(1,:,:);
Ux = Ux(:);
Uy = v(2,:,:);
Uy = Uy(:);
Uz = v(3,:,:);
Uz = Uz(:);
Vu = vecnorm(real(v));
Vu = Vu(:);
tmin = 1;
tmax = size(X,1);
figure(1)
quiver3(X(tmin:tmax),Y(tmin:tmax),Z(tmin:tmax),...
    imag(Ux(tmin:tmax)),imag(Uy(tmin:tmax)),imag(Uz(tmin:tmax)),3)
% hold on
figure(2)
scatter3(X(tmin:tmax),Y(tmin:tmax),Z(tmin:tmax),1e1*(Vu(tmin:tmax)+eps)/max(Vu(tmin:tmax)),Vu(tmin:tmax),'filled')
% figure()
% scatter3(X,Y,Z,10,Ux,'filled')
% figure()
% scatter3(X,Y,Z,10,Uy,'filled')
% figure()
% scatter3(X,Y,Z,10,Uz,'filled')
%tetramesh(TR,'FaceAlpha',0.2)