function feObj = loadFE(FE)
%LOAD Summary of this function goes here
%   Detailed explanation goes here
    s = what('ofem_v2');
    s = [s.path,'/+elements/+preComp/'];
    load([s,FE,'.mat']);
    feObj = fe;
end

