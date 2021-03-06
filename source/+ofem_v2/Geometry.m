classdef Geometry < handle
    %Geometry stores the data of the mesh. (eg. coordinates, parts etc)
    %   Detailed explanation goes here
    
    properties (SetAccess=public)
        % dimension of the space described by the mesh
        dim;

        % ofem.matrixarray of size(co,3) length containing the coordinates
        % of the mesh
        co; 
        
        % Nco is the number of nodes in the loaded mesh. This number is not
        % changed by operations implicitly extending the co
        % ofem.matrixarray.
        Nco;
        
        % ed is a matrix size(ed,2)=2 containing the edge to nodes mapping.
        % More precisely, each row contains the ids of the two nodes that
        % represent the respective edge. The ids are sorted in ascending
        % order.
        ed;
		el2ed;
		Ned;
        
        % refTet switches between one of two reference tetrahedra for the
        % edge-based elements. due to the rotational invariance of edge
        % elements, two reference tetrahedra have to be used
        refTet;
		
		fa;
		el2fa;
		Nfa;
		
		fa2ed;
		Nint;
        
        % el is a matrix with dim+1=size(el,2) containing the
        % element to coordinates mapping. More precisely, each row contains
        % the indices into the co matrix representing the coordinates of
        % its vertices.
        el;
        
        % parts is a cell array with 3=size(parts,1) containing in each
        % column the name of the element set, the index into the material
        % data mat describing the material the element set is build from
        % and a vector of indices into el describing which elements
        % the element set is created from.
        parts;

        % bd is a cell array with 2=size(bd,1) containing in each column
        % the name of the set and a cell array. For dim==2 this cell array
        % is 2x3 and for dim==3 it is 2x4. In any case the first row
        % contains the names of the subset. The second row contains vectors
        % of indices of elements, i.e., indices into el. For dim==2
        % the boundary bd{2,1} is 'The boundary created taking the first
        % edge of the elements with indices bd{2,1}{2,1}, following the 
        % second edge of the elements with indices bd{2,1}{2,2} and the
        % third  edge of the elements with indices bd{2,1}{2,3}'.
        % The procedure is analogous for dim==2, however, a tetrahedron is
        % build up from 4 faces.
        bd;
		
		% roi is a cell array similar to bd
		% here regions of interest are stored
		roi;
        
        type; 

        filetype;
        
        %Jacobian Data
        DinvT;
        detD;
        Dk;
        

    end
    %% 
    
    methods
        
        function obj = Geometry()
        end
        
        function setMaterial (obj, partName, material)
            %cell array to set which matrix needs wich parameter
            if isscalar(partName)
                partIndex = partName;
            else
                [~, number_colums] = size(obj.parts);
                for i= 1:number_colums
                   if isequal(obj.parts{1,i}, partName)
                       partIndex = i;
                       break; 
                   end
                end
            end
            obj.parts{2, partIndex} = material;
        end
        
        function setForce(obj, partName, force)
            if isscalar(partname)
                partIndex = partname;
            else
                [~, number_colums] = size(obj.parts);
                for i= 1:number_colums
                   if isequal(obj.parts{1,i}, partName)
                       partIndex = i;
                       break; 
                   end
                end
            end
            obj.parts{4, partIndex} = force;
        end

      
        function reset(obj)
        %RESET reset the mesh to its initial configuration.
  %TODO is it realy to its inital coonfig?      
                obj.co    = [];
                obj.el    = [];
                obj.bd    = [];
                obj.dim   = [];
                obj.parts = [];
                obj.Nco   = [];
                obj.ed = [];
        end
        
        function hypercube(obj,x0,varargin)
        %hypercube creates a hypercube
        %
        %  hypercube(x0) creates a unit hypercube with lower-left corner at
        %  x0. The dimension of the embedding space is numel(x0).
        %
        %  hypercube(x0,extent) additionally sets the extent of the
        %  hypercube. extent can either be a scalar, in which case the
        %  hypercube has all lengths equal, or must match the dimension of
        %  x0.
        %
        %  In 2D, the mesh consists of two triangles, in 3D it consists of
        %  6 tetrahedra.
        %
        %  hypercube(x0,extent,symmetric) additionally requests to create
        %  a symmetric triangulation, i.e., the created hypercube has an
        %  additional node in the center of gravity. Works for 2D only for
        %  now.
        %
        
            p=inputParser;
            addRequired(p,'x0',@(x) validateattributes(x,{'numeric'},{'nonempty','real','vector'}));
            addOptional(p,'extent',1,@(x) validateattributes(x,{'numeric'},{'nonempty','real','vector'}));
            addOptional(p,'symmetric',false,@(x) validateattributes(x,{'logical'},{'nonempty','scalar'}));
            parse(p,x0,varargin{:});

            extent    = p.Results.extent;
            symmetric = p.Results.symmetric;

            if any(extent<=0)
                error('ofem:mesh:create:argchk','The extent must be positive');
            end

            if numel(extent)==1
                extent=repmat(extent,numel(x0),1);
            else
                if numel(x0)~=numel(extent)
                    error('ofem:mesh:create:argchk','numel(extent) must match numel(x0) or be scalar');
                end
            end

            LL=x0;
            UR=x0+extent;

            obj.reset();

            obj.dim = numel(LL);
            
            switch obj.dim
                case 1
                    % coordinates
                    obj.Nco = 2;
                    obj.co = [LL(1); UR(1)];

                    % elements
                    obj.el = [ 1,2 ];                    

                    % element type
                    obj.type = 'edge';
                    
                    % sidesets
                    obj.bd = {'boundary';{'boundary_P1','boundary_P2'; 1,1}};
                    
                case 2
                    if ~symmetric
                        % coordinates
                        obj.Nco = 4;
                        obj.co  = [         ...
                            [LL(1), LL(2)]; ...
                            [UR(1), LL(2)]; ...
                            [UR(1), UR(2)]; ...
                            [LL(1), UR(2)]  ...
                            ];
                        
                        % elements
                        obj.el = [ ...
                            1,2,4; ...
                            2,3,4  ...
                            ];
                        
                        % element type
                        obj.type = 'tri';
                        
                        % sidesets
                        obj.bd = {'boundary';{'boundary_E1','boundary_E2','boundary_E3'; [1;2],2,1}};
                    else
                        % coordinates
                        MM=0.5*(LL+UR);
                        obj.Nco = 5;
                        obj.co  = [         ...
                            [LL(1), LL(2)]; ...
                            [UR(1), LL(2)]; ...
                            [UR(1), UR(2)]; ...
                            [LL(1), UR(2)]; ...
                            [MM(1), MM(2)]; ...
                            ];

                        % elements
                        obj.el = [ ...
                            1,2,5; ...
                            2,3,5; ...
                            3,4,5; ...
                            4,1,5; ...
                            ];                    

                        % element type
                        obj.type = 'tri';

                        % sidesets
                        obj.bd = {'boundary';{'boundary_E1','boundary_E2','boundary_E3'; [1;2;3;4],[],[]}};
                    end
                    
                case 3
                    % coordinates
                    obj.Nco = 8;
                    obj.co  = [        ...
                        [LL(1), LL(2), LL(3)]; ...
                        [UR(1), LL(2), LL(3)]; ...
                        [UR(1), UR(2), LL(3)]; ...
                        [LL(1), UR(2), LL(3)]; ...
                        [LL(1), LL(2), UR(3)]; ...
                        [UR(1), LL(2), UR(3)]; ...
                        [UR(1), UR(2), UR(3)]; ...
                        [LL(1), UR(2), UR(3)]; ...
                        ];
                    
                    % elements
                    obj.el = [ ...
                        1,2,3,7; ...
                        1,3,4,7; ...
                        1,2,7,6; ...
                        5,1,8,7; ...
                        1,4,8,7; ...
                        5,1,7,6; ...
                        ];
                    
                    % element type
                    obj.type = 'tet';

                    % sidesets
                    % S1: 1, 2, 3
                    % S2: 1, 2, 4
                    % S3: 2, 3, 4
                    % S4: 1, 3, 4
                    obj.bd = {'boundary';{'boundary_S1','boundary_S2','boundary_S3','boundary_S4'; [1;2;4;5],[3;6],[1;2;3;5],[4;6]}};
                    
                otherwise
                    error('ofem:mesh:InvalidMesh',...
                          'hypercube requested for dim>3');
            end

            obj.co = ofem_v2.tools.matrixarray(permute(obj.co,[2,3,1]));

            % parts
            obj.parts = {'hypercube';1;(1:size(obj.el,1))'};

            % materials
            obj.mat{1}.name = 'air';
            obj.mat{1}.ELASTIC = [0;0];
            obj.mat{1}.DENSITY = 0.001;
            obj.mat{1}.CONDUCTIVITY = 0;
            obj.mat{1}.SPECIFIC_HEAT = 1.012;
        end

        function load_from_msh(obj,msh_file_name)
            %LOAD_FROM_MSH loads mesh data from a .msh file
            [filepath,filename,fileext] = fileparts(msh_file_name);
            if fileext==""
                msh_file_name = fullfile(char(filepath),char(strcat(filename,".msh")));
            end
            if filepath==""
                msh_file_name = fullfile(pwd,char(msh_file_name));
            end
            
            % check if msh file exists
            if ~exist(msh_file_name,'file') % for whatever reason this fails on the cluster
                error('ofem:mesh:FileNotFound','msh-file "%s" not found.\n',msh_file_name);
            end
            
            tic
            msh=msh_read_file(msh_file_name);
            obj.filetype = 'msh';
            info.time.file_read=toc;
            
            % msh is a 2x3 cell array
            % {2,1} contains the physicalnames defined in gmsh
            % {2,2} contains the coordinates
            % {2,3} contains the elements
            
            if(isempty(msh{2,3}{2,3}))
                msh{2,2}{2,1}(:,3) = [];
                obj.co = ofem_v2.tools.matrixarray(reshape(msh{2,2}{2,1}',2,1,[]));
                obj.dim = 2;
                obj.type = 'tri';
                obj.el = msh{2,3}{2,2}(:,2:4);
            else
                obj.co = ofem_v2.tools.matrixarray(reshape(msh{2,2}{2,1}',3,1,[]));
                obj.dim = 3;
                obj.type = 'tet';
                obj.el = msh{2,3}{2,3}(:,2:5);
            end
            obj.Nco = size(obj.co,3);
			obj.Nint = size(obj.el,1);
%             obj.el = unique(obj.el(),'rows','legacy');
            
            part = 1;
            bd = 1;
			roi = 1;
            Idx = 1:size(obj.el,1);
            for i=1:size(msh{2,1},2)
                if startsWith(msh{2,1}{1,i},'M:')
                    obj.parts{1,part} = msh{2,1}{1,i}(3:end); % Er schreibt es da nicht rein :'(                 
                    if obj.dim == 2
                        obj.parts{3,part} = Idx(msh{2,3}{2,2}(:,1)==msh{2,1}{2,i});
                        obj.parts{2,part} = 1; %Material extension
                    end
                    if obj.dim == 3
                        obj.parts{3,part} = Idx(msh{2,3}{2,3}(:,1)==msh{2,1}{2,i});
                        obj.parts{2,part} = 1; %Material extension
                    end
                    part = part+1;
                end
                if startsWith(msh{2,1}{1,i},'BD:')
                    obj.bd{1,bd} = msh{2,1}{1,i}(4:end);                   
                    if obj.dim == 2
                        obj.bd{2,bd} = msh{2,3}{2,1}(msh{2,3}{2,1}(:,1)==msh{2,1}{2,i},2:end);
                    end
                    if obj.dim == 3
                        obj.bd{2,bd} = msh{2,3}{2,2}(msh{2,3}{2,2}(:,1)==msh{2,1}{2,i},2:end);
                    end
                    bd = bd+1;
				end
				if startsWith(msh{2,1}{1,i},'ROI0:')
                    obj.roi{1,roi} = msh{2,1}{1,i}(6:end);
                    obj.roi{2,roi} = msh{2,3}{2,1}(msh{2,3}{2,1}(:,1)==msh{2,1}{2,i},2:end);
                    roi = roi+1;
				end
				if startsWith(msh{2,1}{1,i},'ROI1:')
                    obj.roi{1,roi} = msh{2,1}{1,i}(6:end);
                    obj.roi{2,roi} = msh{2,3}{2,1}(msh{2,3}{2,1}(:,1)==msh{2,1}{2,i},2:end);
                    roi = roi+1;
				end
				if startsWith(msh{2,1}{1,i},'ROI2:')
					obj.roi{1,roi} = msh{2,1}{1,i}(6:end);
                    obj.roi{2,roi} = msh{2,3}{2,2}(msh{2,3}{2,2}(:,1)==msh{2,1}{2,i},2:end);
                    roi = roi+1;
                end
                
            end
            jacobiandata(obj);
        end
        
        function export_UCD(obj,folder_name,file_name,meta,varargin)
        %EXPORT_UCD exports solution to inp file.
        %
        % EXPORT_UCD(file_name,meta) exports a solution given by meta to a
        % inp file given by file_name. meta is explained below.
        %
        % EXPORT_UCD(file_name,meta1,meta2,...) exports several solutions.
        %
        % meta is a cell array numel(meta)==3, where meta{2} is the desired
        % function to be exported and must be known for every mesh node.
        % Note that ,e.g., for P2 elements the function must be known on
        % every additional node. meta{2} is not allowed to be sparse.
        % meta{1} is a string containing the name of the function to be
        % exported and meta{3} is a string containinig the name of the unit
        % associated to the function values.
        %
        % The exported inp file uses the UCD format, thus can be loaded,
        % e.g., in ParaView.
        %

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
            Nc   = size(obj.el,1)   ; % number of cells
            Nn   = obj.Nco          ; % number of nodes
%             Nn  = size(obj.co,1);
            Nm   = size(obj.parts,2); % number of materials

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
                                               Nc, ...
                                               size(sol,2),...
                                               size(solC,2));

            % nodes => nodeID <x_1> <x_2> <x_3>
            co_virt=zeros(Nn,3);
            co_virt(:,1:obj.dim)=squeeze(obj.co(:,1,1:Nn))';
            fprintf(fileID,'%d %e %e %e\r\n',[1:Nn;co_virt']);

            % cells => cellID matID type <nodeID_1> ... <nodeID_N>
            currID=0;

            % As of now ParaView, or more precisely VTK, does not support
            % reading higher order elements although, at least, second
            % order elements are supported by UCD.
            switch obj.type
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
        
        function jacobiandata(obj)
        %JACOBIANDATA returns the Jacobians and their determinants
        %
        % [DinvT,detD]=JACOBIANDATA, for a transformation
        % x=\Phi_k(y)=D_ky+n_k
        % from a reference element given in y-coordinates to the global
        % element \Omega_k given in x-coordinates, returns the Jacobians
        % D_k^{-T} and their determinants. DinvT is a matrixarray build up
        % as follows
        %   DinvT(:,:,i) = D_i^{-T}
        % likewise
        %   detD(:,:,i) = det(DinvT(:,:,i))
        %
            switch obj.type
                case 'edge'
                    % edge
                    tmp=obj.co(:,1,obj.el(:,2))-obj.co(:,1,obj.el(:,1));
                    obj.DinvT = 1/tmp;
                    obj.detD  = tmp;

                case 'tri'
                    % triangle
                    e12 = obj.co(:,1,obj.el(:,2))-obj.co(:,1,obj.el(:,1));
                    e13 = obj.co(:,1,obj.el(:,3))-obj.co(:,1,obj.el(:,1));
					
					obj.Dk = [e12,e13];

                    obj.detD = dot(rot(e12),e13,1);

                    obj.DinvT = [ -rot(e13), rot(e12) ]./obj.detD;

                case 'quad'
                    % quadrilateral
                    error('ofem:mesh:InvalidMesh',...
                          'Quadrilateral meshes not supported so far!');

                case 'tet'
                    % tetrahedron
                    e12 = obj.co(:,1,obj.el(:,2))-obj.co(:,1,obj.el(:,1));
                    e13 = obj.co(:,1,obj.el(:,3))-obj.co(:,1,obj.el(:,1));
                    e14 = obj.co(:,1,obj.el(:,4))-obj.co(:,1,obj.el(:,1));
                    
                    obj.Dk = [e12,e13,e14];
                    
                    obj.detD = dot(e12,cross(e13,e14,1),1);

                    obj.DinvT = [ cross(e13,e14,1) , ...
                              cross(e14,e12,1) , ...
                              cross(e12,e13,1) ]./obj.detD;

                case 'hex'
                    % hexahedron
                    error('ofem:mesh:InvalidMesh',...
                          'Hexahedral meshes not supported so far!');

                otherwise
                    error('ofem:mesh:Unspecified',...
                          'Unspecified error found');
            end
		end
		function create_edges(obj)
        %CREATE_EDGES creates edges information and element to edges mapping
        %
            Nt = size(obj.el,1);

            switch obj.type
                case 'edge'
                    %% edge
                    obj.ed        = obj.el;
                    obj.el2ed     = 1:Nt;
                    obj.el2edsign = ones(Nt,1);
                    return;

                case 'tri'
                    %% triangle
                    obj.ed = [obj.el(:,[1 2])'; ...
                              obj.el(:,[1 3])'; ...
                              obj.el(:,[2 3])'];

                case 'quad'
                    %% quadrilateral
                    obj.ed = [obj.el(:,[1 2]); ...
                              obj.el(:,[2 3]); ...
                              obj.el(:,[3 4]); ...
                              obj.el(:,[4 1])];

                case 'tet'
                    %% tetrahedron
                    obj.ed = [obj.el(:,[1 2])'; ...
                              obj.el(:,[1 3])'; ...
                              obj.el(:,[1 4])'; ...
                              obj.el(:,[2 3])'; ...
                              obj.el(:,[2 4])'; ...
                              obj.el(:,[3 4])'];

                case 'hex'
                    %% hexahedron
                    obj.ed = [obj.el(:,[1 2]); ...
                              obj.el(:,[2 3]); ...
                              obj.el(:,[3 4]); ...
                              obj.el(:,[4 1]); ...
                              obj.el(:,[5 6]); ...
                              obj.el(:,[6 7]); ...
                              obj.el(:,[7 8]); ...
                              obj.el(:,[8 5]); ...
                              obj.el(:,[1 5]); ...
                              obj.el(:,[2 6]); ...
                              obj.el(:,[3 7]); ...
                              obj.el(:,[4 8])];

                otherwise
                    error('ofem:mesh:Unspecified',...
                          'Unspecified error found');
            end
            %Nt            = size(obj.ed,1)/2;
            obj.ed        = (reshape(obj.ed,2,[]))';
            obj.ed    = sort(obj.ed,2);
            [obj.ed,~,ic] = unique(obj.ed,'rows','legacy');
            %[obj.ed,J]    = sortrows(obj.ed);
            obj.el2ed     = (reshape(ic,[],Nt))';
            %obj.el2edsign(:,:,:) = 1;
            obj.Ned = size(obj.ed,1);
        end
        
        %%
        function create_faces(obj)
        %CREATE_FACES creates faces information and element to faces mapping
        %
            Nt = size(obj.el,1);

            switch obj.type
                case 'edge'
                    %% edge
                    warning('ofem:mesh:Invalid',...
                            'create_faces called but mesh is composed of edges!');
                    return;

                case 'tri'
                    %% triangle
                    obj.fa        = obj.el;
                    obj.el2fa     = (1:Nt)';

                case 'quad'
                    %% quadrilateral
                    obj.fa        = obj.el;
                    obj.el2fa     = 1:Nt;
%                     obj.el2fasign = ones(Nt,1);

                case 'tet'
                    %% tetrahedron
                    % S1: 1, 2, 3
                    % S2: 1, 2, 4
                    % S3: 2, 3, 4
                    % S4: 1, 3, 4 switched 3 and 4 for nedelec
                    % implementation and DOF handling
                    % Right handed system with outer normal
                    obj.fa = [obj.el(:,[1 2 3]); ...
                              obj.el(:,[1 2 4]); ...
                              obj.el(:,[1 3 4]); ...
                              obj.el(:,[2 3 4]);];

                case 'hex'
                    %% hexahedron
                    error('ofem:mesh:InvalidMesh',...
                          'Hexahedral meshes not supported so far!');

                otherwise
                    error('ofem:mesh:Unspecified',...
                          'Unspecified error found');
			end

            obj.fa    = sort(obj.fa,2);
            [obj.fa,~,ic] = unique(obj.fa,'rows','legacy');
            obj.el2fa     = reshape(ic,Nt,[]);
  			obj.Nfa = size(obj.fa,1);
		end
		
		function connectFa2Ed(obj)
			falist = obj.fa;
			edlist = obj.ed;
			ed1 = falist(:,[1,2]);
			ed2 = falist(:,[2,3]);
			ed3 = falist(:,[1,3]);
			obj.fa2ed = zeros(size(falist));
			[~,obj.fa2ed(:,1)] = ismember(ed1,edlist,'rows');
			[~,obj.fa2ed(:,2)] = ismember(ed2,edlist,'rows');
			[~,obj.fa2ed(:,3)] = ismember(ed3,edlist,'rows');
        end
        
        function reorderAC(obj)
            [~,a] = min(obj.el,[],2);
            a1 = a==1;
            a2 = a==2;
            a3 = a==3;
            a4 = a==4;
            
            idx = [~a1,a2,a3,a4];
            idx(:,1:3) = idx(:,1:3)|circshift(idx(:,1:3),1,2);
            
            while true
                [~,a] = min(obj.el,[],2);
                sLower = a==1;
                
                if sum(~sLower) == 0
                    break
                end
                
                idx(sLower,:) = false;
                
                for i=1:size(obj.el,1)
                    obj.el(i,idx(i,:)) = circshift(obj.el(i,idx(i,:)),1,2);
                end
            end
            
            idx = logical(ones(size(obj.el,1),4));
            idx(:,1) = false;
            
            while true
                [~,b] = max(obj.el,[],2);
                sUpper = b==4;
                if sum(~sUpper) == 0
                    break
                end
                idx(sUpper,:) = false;
                
                for i = 1:size(obj.el,1)
                    obj.el(i,idx(i,:)) = circshift(obj.el(i,idx(i,:)),1,2);
                end
            end
            
            obj.refTet = double(obj.el(:,2)>obj.el(:,3))+1;
            
        end
     
    end
end








