function getFieldData(res,roi,file)
%GETFIELDDATA induced field on a region of interest surface (e.g. as
%input for inverseFEM)
%
%   getFieldData(res)                       first ROI, fieldData.mat
%   getFieldData(res, 'Upper', 'fieldData.mat')
%
%   Saves the elements adjacent to the ROI faces, their centres x_given
%   and the induced field b_given.
    if nargin < 2 || isempty(roi), roi = 1; end
    if nargin < 3, file = 'fieldData.mat'; end
    mesh = res.mesh;
    if ischar(roi) || isstring(roi)
        roi = find(strcmp(mesh.roi(1,:),roi),1);
    end
    [~,faces] = ismember(sort(mesh.roi{2,roi},2),mesh.fa,'rows');
    [~,elements] = ismember(faces,mesh.el2fa);
    b_given = res.Bind(:,elements); %#ok<NASGU>
    x_given = res.X(:,elements); %#ok<NASGU>
    save(file,'elements','x_given','b_given');
end
