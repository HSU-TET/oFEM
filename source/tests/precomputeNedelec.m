s = what('ofem_v2');
for i = 0:4
    fe = ofem_v2.elements.NedelecZagl(3,i);
    [N,curlN] = fe.computeBasis;
    save([s.path,sprintf('/+elements/+preComp/NE%d.mat',i)],'fe');
end