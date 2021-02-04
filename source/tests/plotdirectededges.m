v1 = mesh.co(:,1,mesh.el(:,2))-mesh.co(:,1,mesh.el(:,1));
v2 = mesh.co(:,1,mesh.el(:,3))-mesh.co(:,1,mesh.el(:,1));
v3 = mesh.co(:,1,mesh.el(:,4))-mesh.co(:,1,mesh.el(:,1));
v4 = mesh.co(:,1,mesh.el(:,3))-mesh.co(:,1,mesh.el(:,2));
v5 = mesh.co(:,1,mesh.el(:,4))-mesh.co(:,1,mesh.el(:,2));
v6 = mesh.co(:,1,mesh.el(:,4))-mesh.co(:,1,mesh.el(:,3));

v = [squeeze(v1),squeeze(v2),squeeze(v3),squeeze(v4),squeeze(v5),squeeze(v6)];

start1 = mesh.co(:,1,mesh.el(:,1));
start2 = mesh.co(:,1,mesh.el(:,1));
start3 = mesh.co(:,1,mesh.el(:,1));
start4 = mesh.co(:,1,mesh.el(:,2));
start5 = mesh.co(:,1,mesh.el(:,3));
start6 = mesh.co(:,1,mesh.el(:,4));

start = [squeeze(start1),squeeze(start2),squeeze(start3),squeeze(start4),squeeze(start5),squeeze(start6)];

%quiver3(start(1,:),start(2,:),start(3,:),v(:,1),v(:,2),v(:,3))
%plot3(start(1,:)+v(1,:)/2,start(2,:)+v(2,:)/2,start(3,:)+v(3,:)/2,'x')

for i = 1:size(start,2)
    quiver3(start(1,i),start(2,i),start(3,i),v(1,i),v(2,i),v(3,i),0.5);
    hold on
end