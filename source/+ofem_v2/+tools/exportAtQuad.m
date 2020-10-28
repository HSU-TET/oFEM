function exportAtQuad(folder_name,file_name,X,meta,varargin)
%EXPORTATQUAD Summary of this function goes here
%   Detailed explanation goes here
                % create output file________________________________________________
            if ~exist(folder_name,'dir')
                warning('ofem:mesh:FolderNotExistent',...
                        'The given folder is not existent, still I''ll create it.');
                mkdir(folder_name);
            end

            [pathstr,file_name,ext]=fileparts(file_name);
            if ~isempty(pathstr) % user requested current working directory
                warning('ofem:mesh:InvalidArgument',...
                        'file_name contains a path. I''ll ignore it');
            end
            if isempty(ext)
                file_name=strcat(file_name,'.inp');
            end
            file_name=fullfile(folder_name,file_name);


            % basic computations________________________________________________
            Nf   = nargin-3         ; % number of functions
%            Nc   = size(obj.el,1)   ; % number of cells
            Nn   = size(X,1)          ; % number of nodes
%             Nn  = size(obj.co,1);
%            Nm   = size(obj.parts,2); % number of materials
            type = 'tet';

            % 
            name{1} = 0;
            unit{1} = 0;
            Nfd  = 0;
            nameC{1} = 0;
            unitC{1} = 0;
            NfdC  = 0;
            solC = [];
            sol = [];
            j = 0;
            k = 0;

            if issparse(meta{2})
                warning('ofem:mesh:InvalidArgument',...
                        'u is not allowed to be sparse. I''ll convert it to full');
            end
            if numel(meta)~=4&&numel(meta)~=3
                error('ofem:mesh:InvalidArgument',...
                      'meta is expected to be a cell array with numel(meta)==4');
            end
            if numel(meta)==3
                name{1} = meta{1};
                sol     = full(meta{2});
                Nfd (1) = size(meta{2},2);
                unit{1} = meta{3};
                j = j+1;
            elseif meta{4} =='Cell'
                nameC{1} = meta{1};
                solC     = full(meta{2});
                NfdC(1)  = size(meta{2},2);
                unitC{1} = meta{3};
                k = k+1;
            end
            

            for i=1:nargin-4
                if issparse(varargin{i}{2})
                error('ofem:mesh:InvalidArgument',...
                      'u is not allowed to be sparse. I''ll convert it to full');
                end
                if numel(varargin{i})~=4&&numel(varargin{i})~=3
                    error('ofem:mesh:InvalidArgument',...
                          'meta is expected to be a cell array with numel(meta)==4');
                end

                if numel(varargin{i})==3
                    name{j+1} = varargin{i}{1};
                    Nfd (j+1) = size(varargin{i}{2},2);
                    unit{j+1} = varargin{i}{3};
                    sol(:,end+(1:Nfd(j+1))) = full(varargin{i}{2});
                    j = j+1;
                elseif varargin{i}{4} == 'Cell'
                    nameC{k+1} = varargin{i}{1};
                    NfdC(k+1)  = size(varargin{i}{2},2);
                    unitC{k+1} = varargin{i}{3};
                    solC(:,end+(1:NfdC(k+1))) = full(varargin{i}{2});
                    k = k+1;
                end

            end


            % begin of export________________________________________________
            tic
            fprintf('Starting data export ... ');
            fileID = fopen(file_name,'w+');

            % UCD format

            % header => #nodes #cells #nodedata #celldata #classdata
            fprintf(fileID,'%d %d %d %d 0 \r\n',Nn, ...
                                               size(sol,2),...
                                               size(solC,2));

            % nodes => nodeID <x_1> <x_2> <x_3>
            %co_virt=zeros(X,3);
            co_virt = X;
            fprintf(fileID,'%d %e %e %e\r\n',[1:Nn;co_virt']);

            % cells => cellID matID type <nodeID_1> ... <nodeID_N>
            currID=0;

            % As of now ParaView, or more precisely VTK, does not support
            % reading higher order elements although, at least, second
            % order elements are supported by UCD.
            switch type
                case 'line'
                    formatStr = '%d %d line %d %d\r\n';
                    c2nidx = 1:2;
                case 'tri'
                    formatStr = '%d %d tri %d %d %d\r\n';
                    c2nidx = 1:3;
                case 'tet'
                    formatStr = '%d %d tet %d %d %d %d\r\n';
                    c2nidx = 1:4;
                case 'quad'
                    formatStr = '%d %d quad %d %d %d %d\r\n';
                    c2nidx = 1:4;
                case 'hex'
                    formatStr = '%d %d hex %d %d %d %d %d %d %d %d\r\n';
                    c2nidx = 1:8;
                otherwise
                    error('ofem:mesh:InvalidMesh',...
                          'Unsupported mesh found. Type is %s',obj.type);
            end
        
%             formatStr = ['%d %d ',obj.type,repmat(' %d',1,size(obj.el,2)),'\r\n'];
            
            for i=1:Nm
                partIDs  = obj.parts{3,i};
				partIDs = reshape(partIDs,1,[]);
                NpartIDs = numel(partIDs);
                fprintf(fileID,formatStr,[partIDs; i*ones(1,NpartIDs); obj.el(partIDs,c2nidx)']);
                currID=currID+NpartIDs;
            end

            % write functions________________________________________________
            % header => #func #colsfunc_1 ... #colsfunc_N


            % values => nodeID funcval_1 ... funcval_N
            if Nfd~=0
                fprintf(fileID,'%d ',[size(Nfd,2), Nfd]);
                fprintf(fileID,'\r\n');

                % names, units
                for i=1:size(Nfd,2)
                    fprintf(fileID,'%s, %s\r\n',name{i},unit{i});
                end
                formatStr = ['%d ',repmat(' %e',1,size(sol,2)),'\r\n'];
                fprintf(fileID,formatStr,[1:Nn;sol(1:Nn,:)']);
            end
            
            if NfdC ~=0
                fprintf(fileID,'%d ',[size(NfdC,2), NfdC]);
                fprintf(fileID,'\r\n');

                % names, units
                for i=1:size(NfdC,2)
                    fprintf(fileID,'%s, %s\r\n',nameC{i},unitC{i});
                end
                formatStr = ['%d ',repmat(' %e',1,size(solC,2)),'\r\n'];
                fprintf(fileID,formatStr,[1:Nc;solC(1:Nc,:)']);
            end

            % close file print info and exit________________________________________________
            fclose(fileID);

            t=toc;
            fprintf('done. t=%f s\n',t);
    
    
end

