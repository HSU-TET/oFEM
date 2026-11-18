function exportResults(res,name,folder)
%EXPORTRESULTS UCD export (cell data) of a movingConductor result
%
%   exportResults(res, name)              into ./export
%   exportResults(res, name, folder)
    if nargin < 3, folder = fullfile(pwd,'export'); end
    res.mesh.export_UCD(folder,name, ...
        {'B_tot',res.B','','Cell'},{'A',res.A','','Cell'}, ...
        {'B_ind',res.Bind','','Cell'},{'J',res.J','','Cell'}, ...
        {'J_ind',res.Jind','','Cell'});
end
