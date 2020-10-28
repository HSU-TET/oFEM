k = zeros(3,6);
for i = 0:2
    tic
    k(:,i+1) = fichera(i);
    time(i+1) = toc
end
sum(time)