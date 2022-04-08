%% Use this functions to recompute the precomputed Basis functions

s = what('ofem_v2');
s = [s.path,'/+elements/+preComp/'];

% H1
for dim = 1:3
	for deg = 1:3
		fe = ofem_v2.elements.H1Element(dim,deg);
		fe.computeBasis;
		save([s,'H1_',num2str(dim),'D_','Order_',num2str(deg),'.mat'],'fe');
	end
end

% H Curl
for dim = 2:3
	for deg = 0:3
		fe = ofem_v2.elements.HCurlElement(dim,deg);
		fe.computeBasis;
		save([s,'HCurl_',num2str(dim),'D_','Order_',num2str(deg),'.mat'],'fe');
	end
end