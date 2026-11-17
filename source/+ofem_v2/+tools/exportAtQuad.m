function exportAtQuad(folder_name,file_name,X,meta,varargin)
%EXPORTATQUAD exports values at arbitrary points (e.g. quadrature points)
%as a UCD point cloud (.inp, readable by ParaView)
%
%   exportAtQuad(folder_name, file_name, X, {name, values, unit}, ...)
%
%   X       points, N x dim, or dim x Ne x Np as returned by
%           ofem_v2.tools.reconstruct
%   values  N x k, or k x Ne x Np (same layout as X)
%
%   Every point becomes a 'pt' cell, the values are written as node data.
%   (The previous version referenced undefined variables and wrote an
%   invalid header, so it could not run.)

    if ~exist(folder_name,'dir')
        mkdir(folder_name);
    end
    [pathstr,file_name,ext] = fileparts(file_name);
    if ~isempty(pathstr)
        warning('ofem:mesh:InvalidArgument','file_name contains a path. I''ll ignore it');
    end
    if isempty(ext)
        ext = '.inp';
    end
    file_name = fullfile(folder_name,[file_name,ext]);

    P = toRows(X);
    Nn = size(P,1);
    co = zeros(Nn,3);
    co(:,1:size(P,2)) = P;

    fields = [{meta},varargin];
    names = cell(1,numel(fields));
    units = names;
    cols = zeros(1,numel(fields));
    sol = zeros(Nn,0);
    for i = 1:numel(fields)
        f = fields{i};
        if numel(f) < 3
            error('ofem:exportAtQuad:InvalidArgument','meta must be {name, values, unit}.');
        end
        v = full(toRows(f{2}));
        if size(v,1) ~= Nn
            error('ofem:exportAtQuad:InvalidArgument', ...
                'Field %s has %d values for %d points.',f{1},size(v,1),Nn);
        end
        names{i} = f{1};
        units{i} = f{3};
        cols(i) = size(v,2);
        sol = [sol, v]; %#ok<AGROW>
    end

    fid = fopen(file_name,'w');
    fprintf(fid,'%d %d %d 0 0\r\n',Nn,Nn,size(sol,2));
    fprintf(fid,'%d %e %e %e\r\n',[1:Nn; co']);
    fprintf(fid,'%d 1 pt %d\r\n',[1:Nn; 1:Nn]);
    if ~isempty(sol)
        fprintf(fid,'%d',numel(fields));
        fprintf(fid,' %d',cols);
        fprintf(fid,'\r\n');
        for i = 1:numel(fields)
            fprintf(fid,'%s, %s\r\n',names{i},units{i});
        end
        fprintf(fid,['%d',repmat(' %e',1,size(sol,2)),'\r\n'],[1:Nn; sol']);
    end
    fclose(fid);
end

function A = toRows(A)
    % N x k from N x k, or k x Ne x Np (k <= 3) to (Ne*Np) x k
    A = double(A);
    if ndims(A) == 3
        A = reshape(A,size(A,1),[])';
    end
end
