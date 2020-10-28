k = [];
for i = [2:4]
    k = [k,thickL(i)];
end
%%
ref = [9.63972384472;11.3452262252;13.4036357679;15.1972519265;19.5093282458;
           19.7392088022;19.7392088022;19.7392088022];
res = k;
res(res<1) = 0;
res(isnan(res)) = 0;
res = sort(res,1);
res(res==0) = [];
res = [res,ref];
%res = abs(res-ref)./ref;
plot(res')
legend()