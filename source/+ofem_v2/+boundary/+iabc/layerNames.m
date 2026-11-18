function layers = layerNames(mesh)
%LAYERNAMES names of the IABC shells 'Layer1', 'Layer2', ... in the mesh,
%sorted from the inside out
    names = mesh.parts(1,:);
    tok = regexp(names,'^Layer(\d+)$','tokens','once');
    has = ~cellfun(@isempty,tok);
    num = cellfun(@(t) str2double(t{1}),tok(has));
    [~,o] = sort(num);
    layers = names(has);
    layers = layers(o);
    if isempty(layers)
        error('iabc:layers','The mesh has no parts ''Layer1'', ''Layer2'', ...');
    end
end
