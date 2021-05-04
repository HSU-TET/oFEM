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
        
        % mat is a cell array containing the materials the meshed model is
        % build up from. Each entry in the array is a structure with a
        % field name containing the name of the material. Dependent on the
        % material there can be additional fields, all uppercase with
        % whitespaces replaced by underscores, containing the desired
        % property as, e.g., DENSITY. parts{2,:} containg the indices
        % into mat.
        mat;
        

    end
   
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
        
        function uniform_refine(obj)
        %uniform_refine uniformly refines the mesh
        %
            if isempty(obj.type) == true
                error('ofem:mesh:InvalidMesh','ofem.mesh needs to hold a valid triangulation');
            end

%             recover_fe=0;
% 
%             if ~isempty(obj.fe)
%                 fe_tmp=obj.fe;
%                 obj.reset();
%                 recover_fe=1;
%             end
            

            [obj.co,elem] = obj.create_midpoints();

            obj.Nco = size(obj.co,3);

            Nel = size(obj.el,1);

            switch obj.type
                case 'edge'
                    %% edge
                    % update elements
                    obj.el = [         ...
                        elem(:,[1,3]); ...
                        elem(:,[3,2])  ...
                        ];

                    % update parts
                    obj.parts(3,:)=cellfun(@(x) [x;x+Nel],obj.parts(3,:),'UniformOutput',false);

                    % update boundaries
                    for i=1:numel(obj.bd(2,:))
                        obj.bd{2,i}{2,2}=obj.bd{2,i}{2,2}+Nel;
                    end

                case 'tri'
                    %% triangle
                    % update elements
                    obj.el = [ ...
                        elem(:,[1,4,6]); ...
                        elem(:,[4,2,5]); ...
                        elem(:,[6,5,3]); ...
                        elem(:,[4,5,6])  ...
                        ];

                    % update parts
                    obj.parts(3,:)=cellfun(@(x) [x;x+Nel;x+2*Nel;x+3*Nel],obj.parts(3,:),'UniformOutput',false);

                    % update boundaries
                    for i=1:numel(obj.bd(2,:))
                        obj.bd{2,i}{2,1}=[obj.bd{2,i}{2,1}      ; obj.bd{2,i}{2,1}+  Nel];
                        obj.bd{2,i}{2,2}=[obj.bd{2,i}{2,2}+  Nel; obj.bd{2,i}{2,2}+2*Nel];
                        obj.bd{2,i}{2,3}=[obj.bd{2,i}{2,3}+2*Nel; obj.bd{2,i}{2,3}      ];
                    end

                case 'tet'
                    %% tetrahedron
%                     error('ofem:mesh:NotImplemented',...
%                           'Tetrahedral meshes not supported so far!');

                    obj.el = [ ...
                        elem(:,[1,5, 6, 7]); ...
                        elem(:,[5,2, 8, 9]); ...
                        elem(:,[6,8, 3,10]); ...
                        elem(:,[7,9,10, 4]); ...
                        elem(:,[7,5, 6, 9]); ...
                        elem(:,[5,8, 6,10]); ...
                        elem(:,[6,8,10, 9]); ...
                        elem(:,[7,9, 6,10]); ...
                        ];

                    % update parts
                    obj.parts(3,:)=cellfun(@(x) [x;x+Nel;x+2*Nel;x+3*Nel;x+4*Nel;x+5*Nel;x+6*Nel;x+7*Nel],obj.parts(3,:),'UniformOutput',false);

                    % update boundaries
                    % S1: 1, 2, 3
                    % S2: 1, 2, 4
                    % S3: 2, 3, 4
                    % S4: 1, 3, 4
                    for i=1:numel(obj.bd(2,:))
                        obj.bd{2,i}{2,1}=[obj.bd{2,i}{2,1}      ; obj.bd{2,i}{2,1}+1*Nel; obj.bd{2,i}{2,1}+2*Nel; obj.bd{2,i}{2,1}+5*Nel];
                        obj.bd{2,i}{2,2}=[obj.bd{2,i}{2,2}      ; obj.bd{2,i}{2,2}+1*Nel; obj.bd{2,i}{2,2}+3*Nel; obj.bd{2,i}{2,2}+4*Nel];
                        obj.bd{2,i}{2,3}=[obj.bd{2,i}{2,3}+1*Nel; obj.bd{2,i}{2,3}+2*Nel; obj.bd{2,i}{2,3}+3*Nel; obj.bd{2,i}{2,3}+6*Nel];
                        obj.bd{2,i}{2,4}=[obj.bd{2,i}{2,4}      ; obj.bd{2,i}{2,4}+2*Nel; obj.bd{2,i}{2,4}+3*Nel; obj.bd{2,i}{2,4}+7*Nel];
                    end

                otherwise
                    error('ofem:mesh:Unspecified',...
                          'Unspecified error found');
            end

%             if recover_fe==1
%                 obj.assign_fe(fe_tmp);
%             end
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
        
        function info=load_from_inp(obj,inp_file_name)
        % LOAD_FROM_INP(inp_file_name) loads the inp-file specified by
        % inp_file_name into the ofem_v2.Geometry class.
        %

            % check if extension was specified
            [filepath,filename,fileext] = fileparts(inp_file_name);
            if isempty(fileext)
                inp_file_name = fullfile(filepath,strcat(filename,'.inp'));
            end
            if isempty(filepath)
                inp_file_name = fullfile(pwd,inp_file_name);
            end

            % check if inp file exists
            if ~exist(inp_file_name,'file') % for whatever reason this fails on the cluster
                error('ofem:mesh:FileNotFound','inp-file "%s" not found.\n',inp_file_name);
            end

            tic
            inp=inp_read_file(inp_file_name);
            obj.filetype = 'inp';
            info.time.file_read=toc;

            % inp is a cell array of the following structure
            % inp{:,1}: coordinates
            % inp{:,2}: elements
            % inp{:,3}: nodesets
            % inp{:,4}: sidesets
            % inp{:,5}: properties
            % inp{:,6}: materials

            tic

            %% coordinates
            % size(obj.co)=[Nd,1,Nco], always column vectors
            obj.co  = ofem_v2.tools.matrixarray(permute(inp{2,1}{2,1},[2,3,1]));
            obj.Nco = size(obj.co,3);
            obj.dim = size(obj.co,1);

            %% elements
            elaux  = inp{2,2};
            %% parts
            Nparts    = numel(elaux(1,:));
            Nelems =0;
            % ids       =cell(1,Nparts);
            obj.parts = cell(3,Nparts);
            obj.parts(1,:) = elaux(1,:);
            
            
            for i =1:Nparts
               % ids            = elaux{2,i}(:,1);
                obj.el(elaux{2,i}(:,1),:)  = elaux{2,i}(:,2:end);
                obj.parts{3,i} = elaux{2,i}(:,1);
                Nelems         = Nelems + size(elaux{2,i},1);
            end
            
%             
%             obj.parts = cell(3,Nparts);
% 
%             set name
%             obj.parts(1,:) = elaux(1,:);
%             and indices
%             max_idx=0;
%             for i=1:Nparts
%                 Nel            = size(elaux{2,i},1);
%                 obj.parts{3,i} = max_idx+(1:Nel)';
%                 max_idx        = max_idx+Nel;
%             end
            clear elaux; % ids;

            %% element type
            nodes_per_elem=size(obj.el,2);
            switch obj.dim
                case 1
                    if nodes_per_elem~=2
                        error('ofem:mesh:InvalidMesh',...
                              'In 1D the elements are expected to contain two nodes, but %d nodes were found.',...
                              nodes_per_elem);
                    end
                    obj.type = 'line';

                case 2
                    switch nodes_per_elem
                        case 3
                            obj.type = 'tri';
                        case 4
                            obj.type = 'quad';
                        otherwise
                            error('ofem:mesh:InvalidMesh',...
                              'In 2D the elements are expected to contain three of four nodes, but %d nodes were found.',...
                              nodes_per_elem);
                    end

                case 3
                    switch nodes_per_elem
                        case 4
                            obj.type = 'tet';
                        case 8
                            obj.type = 'hex';
                        otherwise
                            error('ofem:mesh:InvalidMesh',...
                              'In 3D the elements are expected to contain four of eight nodes, but %d nodes were found.',...
                              nodes_per_elem);
                    end
            end

            % post-processing
            obj.jacobiandata;
            obj.create_edges;
            obj.create_faces;
            obj.connectFa2Ed;
            
            %% nodesets
%             bd_ns_name=inp{2,3};

            %% sidesets
            obj.bd=inp{2,4};
           
            % compute the edges for the specified boundary
            for idx = 1:size(inp{2,4},2)
                ss    = obj.bd(2,idx);
                Nss   = size(ss,2);
                eID   = cell(Nss,1);
                switch obj.filetype
                    case 'inp'
                        switch obj.type
                            
                            case 'tri'
                                for i=1:Nss
%                                         eID{i}=unique([obj.el(ss{1,i}{2,1},[1,2]); ...
%                                                        obj.el(ss{1,i}{2,2},[2,3]);...
%                                                        obj.el(ss{1,i}{2,3},[3,1])],'stable');

                                            eID{i}=[obj.el(ss{1,i}{2,1},[1,2]); ...
                                                       obj.el(ss{1,i}{2,2},[2,3]);...
                                                       obj.el(ss{1,i}{2,3},[3,1])];
                                end
                                
                                obj.bd{2,idx} = eID{:};
                        
                            case 'tet'
                                for i=1:Nss
%                                     eID{i} = unique([obj.el2ed(ss{1,i}{2,1},[1,2,4]);...
%                                                   obj.el2ed(ss{1,i}{2,2},[1,3,5]);...
%                                                   obj.el2ed(ss{1,i}{2,3},[4,5,6]);...
%                                                   obj.el2ed(ss{1,i}{2,4},[2,3,6])]);

                                    eID{i} = [obj.el2ed(ss{1,i}{2,1},[1,2,4]);...
                                                  obj.el2ed(ss{1,i}{2,2},[1,3,5]);...
                                                  obj.el2ed(ss{1,i}{2,3},[4,5,6]);...
                                                  obj.el2ed(ss{1,i}{2,4},[2,3,6])];
                                              
                                    


                                end
                                
                                [~,fID] = ismember(sort(eID{:},2),sort(obj.fa2ed,2),'rows');
                                        obj.bd{2,idx} = obj.fa(fID,:);
                            otherwise
                                error('ofem:mesh:dirichletEdges',...
                                      'Only Tetrahedron implemented so far!')
                        end
                    case 'msh'
                        ss    = ss{1};
                        switch obj.type
                            case 'tet'
                                %ss = unique(sort(ss,2),'rows','legacy');
                                ss = [ss(:,2),ss(:,1);ss(:,3),ss(:,2);ss(:,1),ss(:,3)];
                                [~,idx] = ismember(ss,obj.ed,'rows');
                                idx(idx==0) = [];
                                eID{1} = idx;

                            otherwise
                                error('ofem:mesh:dirichletEdges',...
                                      'Only Tetrahedron implemented so far!')
                        end
                    otherwise
                        error('ofem:mesh:dirichletEdges',...
                              'Something went horribly wrong!');
                end
       
%                 obj.bd{2,idx} = eID{:};
                
            end
            

            %% materials
            mataux=inp{2,6};
            if isempty(mataux)
              % if no material present => set air as default material
              obj.mat{1}.name = 'air';
              obj.mat{1}.ELASTIC = [0;0];
              obj.mat{1}.DENSITY = 0.001;
              obj.mat{1}.CONDUCTIVITY = 0;
              obj.mat{1}.SPECIFIC_HEAT = 1.012;

              for i=1:Nparts
                obj.parts{2,i}=1;
              end
            else
              % else read materials
              Nmat=size(mataux,2);
              obj.mat=cell(1,Nmat);
              
              for i=1:Nmat
                obj.mat{i}.name=mataux{1,i};
                Nmatpar=numel(mataux{2,i}(1,:));
                for j=1:Nmatpar
                  obj.mat{i}.(strrep(mataux{2,i}{1,j},' ','_'))=mataux{2,i}{2,j};
                end
                
              end
              clear mataux;
              
              % and associate them with the element sets
              props=inp{2,5};
              Nprops=numel(props(1,:));
              
              matnames = cellfun(@(c) c.name,obj.mat,'UniformOutput',0);
              elnames  = obj.parts(1,:);
              
              for i=1:Nprops
                elidx  = strcmp(props{2,i}{1},elnames );
                matidx = find(strcmp(props{2,i}{2},matnames));
                
                obj.parts{2,elidx}=matidx;
              end
            end
            
            info.time.post_proccess=toc;

            obj.Nint = size(obj.el,1);
            
            clear inp;
            
            
            
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
        
        function [co,el]=create_midpoints(obj)
            Nel=numel(obj.el(:,1));
            
            co = obj.co;

            switch obj.type
                case 'edge'
                    %% edge
                    % create new nodes
                    co2 = (obj.co(:,:,obj.el(:,1))+obj.co(:,:,obj.el(:,2)))/2;
                    co2 = permute(co2,[3,1,2]); % coordinate per row
                    
                    % get rid of doubles
                    [co2,~,ib]=unique(co2,'rows','stable');
                    
                    % add new indices
                    Nco2    = size(co2,1);
                    idx2    = obj.Nco+(1:Nco2);
                    
                    % extend nodes
                    co(:,:,end+1:end+Nco2) = permute(co2,[2,3,1]);
                    
                    % append new nodes to elements
                    el=[obj.el, reshape(idx2(ib),Nel,[])];
                    
                case 'tri'
                    %% triangle
                    % create new nodes
                    mpt1 = permute(obj.co(:,:,obj.el(:,1))+obj.co(:,:,obj.el(:,2)),[3,1,2]);
                    mpt2 = permute(obj.co(:,:,obj.el(:,2))+obj.co(:,:,obj.el(:,3)),[3,1,2]);
                    mpt3 = permute(obj.co(:,:,obj.el(:,3))+obj.co(:,:,obj.el(:,1)),[3,1,2]);
                    
                    co2 = [ mpt1;  mpt2;  mpt3 ]/2;
                    
                    % get rid of doubles
                    [co2,~,ib]=unique(co2,'rows','stable');
                    
                    % add new indices
                    Nco2    = size(co2,1);
                    idx2=obj.Nco+(1:Nco2);
                    
                    % extend nodes
                    co(:,:,end+1:end+Nco2) = permute(co2,[2,3,1]);
                    
                    % append new nodes to elements
                    el=[obj.el, reshape(idx2(ib),Nel,[])];
                    
                case 'quad'
                    %% quadrilateral
                    error('ofem:mesh:NotImplemented',...
                          'Quadrilateral meshes not supported so far!');
                    
                case 'tet'
                    %% tetrahedron
                    % create new nodes
                    mpt1 = permute(obj.co(:,:,obj.el(:,1))+obj.co(:,:,obj.el(:,2)),[3,1,2]);
                    mpt2 = permute(obj.co(:,:,obj.el(:,1))+obj.co(:,:,obj.el(:,3)),[3,1,2]);
                    mpt3 = permute(obj.co(:,:,obj.el(:,1))+obj.co(:,:,obj.el(:,4)),[3,1,2]);
                    mpt4 = permute(obj.co(:,:,obj.el(:,2))+obj.co(:,:,obj.el(:,3)),[3,1,2]);
                    mpt5 = permute(obj.co(:,:,obj.el(:,2))+obj.co(:,:,obj.el(:,4)),[3,1,2]);
                    mpt6 = permute(obj.co(:,:,obj.el(:,3))+obj.co(:,:,obj.el(:,4)),[3,1,2]);
                    
                    co2 = [ mpt1; mpt2; mpt3; mpt4; mpt5; mpt6 ]/2;
                    
                    % get rid of doubles
                    [co2,~,ib]=unique(co2,'rows','stable');
                    
                    % add new indices
                    Nco2    = size(co2,1);
                    idx2=obj.Nco+(1:Nco2);
                    
                    % extend nodes
                    co(:,:,end+1:end+Nco2) = permute(co2,[2,3,1]);
                    
                    % append new nodes to elements
                    el=[obj.el, reshape(idx2(ib),Nel,[])];
                    
                case 'hex'
                    %% hexahedron
                    error('ofem:mesh:NotImplemented',...
                          'Hexahedral meshes not supported so far!');
                    
                otherwise
                    error('ofem:mesh:Unspecified',...
                          'Unspecified error found');
            end
        end
        
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
        
        function [idx, tr, bary] = point_location(obj, xq, opt)
        %point_location returns the index of the element which contains the
        %query-point. idx is the index, tr is the computed aabb-tree and bary 
        %are the baryzentric coordinates of the querypoint, depending on the 
        %vertices of the element in which xq lays. If there are multiple 
        %query-points, the vector xq needs the shape NxM where N is the 
        %number of query-points and M is the dimension of the query-points. 
        %2D and 3D are supported. With varargin the user can give an already
        %computed aabb-tree or he can set the parameter which controls the 
        %dividing of the global boundingboxes. 
        %If the function returns NaN, it means that the query-point is not
        %in the given mesh. 
           
           
		   if isfield(opt,'tr')
			   tr = opt.tr;
		   else
			   tr = [];
		   end
		   if isfield(opt,'op')
			   op = opt.op;
		   else
			   op = [];
		   end
		   tt = [];
		   elidx = [];
		   if isfield(opt,'parts')
			   for i=1:length(opt.parts)
					tt  = [tt;obj.el(obj.parts{3,opt.parts(i)},:)];
					elidx = [elidx,obj.parts{3,opt.parts(i)}];
			   end
		   else
			   tt = obj.el;
		   end
		   pp  = double(reshape(permute(obj.co,[3,1,2]),[],size(obj.co,1)));
           tol = 1e-8;
           tp  = []; 
           tj  = []; 
           %tr  = []; 
           %op  = [];               
		   
           %if (nargin >= +5), tr = varargin{1}; end
           %if (nargin >= +6), op = varargin{2}; end
                             
           if (isempty(tr))
                %compute aabb's for triangles
                bi = pp(tt(:,1),:); bj = pp(tt(:,1),:);
                for ii = 2 : size(tt,2)    
                    bi = min(bi,pp(tt(:,ii),:)) ;
                    bj = max(bj,pp(tt(:,ii),:)) ;
                end
 
                bb = [bi-tol,bj+tol];
                %compute aabb-tree for aabb's
                if (isempty(op))                % scale tree with |pj|
                    op.nobj = ceil(size(tt,1)/...
                                   size(xq,1)) * +4 ;                   
                    op.nobj = max( +32,op.nobj);  
                    op.nobj = min(+256,op.nobj);
                end
                
                tr = maketree(bb,op);           % compute aabb-tree
                
           else
            %check existing aabb-tree
                if (~isfield(tr,'xx') || ~isfield(tr,'ii') || ~isfield(tr,'ll') )
                      error('findtria:incorrectInputClass','Incorrect aabb-tree.') ;
                end
           end
                
           [tm,~] = scantree(tr,xq,@partpts);  
           
           function [j1,j2] = partpts (pi,b1,b2)
                %PARTPTS partition points between boxes B1,B2 for SCANTREE.

                j1 = true(size(pi,1),1);
                j2 = true(size(pi,1),1);

                nd = size(b1,2) / +2;

                for ax = +1 : nd
                %remains TRUE if inside bounds along axis AX
                j1 = j1 & pi(:,ax)>=b1(ax+nd*0) ...
                        & pi(:,ax)<=b1(ax+nd*1) ;

                %remains TRUE if inside bounds along axis AX
                j2 = j2 & pi(:,ax)>=b2(ax+nd*0) ...
                        & pi(:,ax)<=b2(ax+nd*1) ;    
                end

            end
           
           
           
           ic     = cell(size(tm.ii,1),1);
           jc     = cell(size(tm.ii,1),1);
           for ip = 1 : size(tm.ii,1)
           %extract trias/items per tile and search 
                ni      = tm.ii(ip,1) ;          % node (in tree)
                [pi,ti] = testpts(pp,tt,xq,tm.ll{ip,1},tr.ll{ni,1});
                ic{ip}  = pi;
                jc{ip}  = ti;        
           end 
                
           %concat matches into arrays
           ii = vertcat(ic{:});
           tj = vertcat(jc{:});
           nj = length(tj);
           %form sparse-style indexing
           [ii,ix] = sort (ii) ; 
           tj      = tj(ix);
           ix      = find(diff(ii) > +0);
           %IN = TRUE if we found any matches
           in     = false(size(xq,1),1);
           in(ii) = true;
           %ptrs into TJ for each point in PJ
           tp       = zeros(size(xq,1),2);
           tp(in,1) = [+1; ix+1];
           tp(in,2) = [ix; nj+0];

           tr.bb = bb;  

           idx     = nan(size(tp,1),1);
           in      = tp(:,1) > +0;
           idx(in) = tj(tp(in,+1));
        
        
           function [tr] = maketree(rp,varargin) 
            
                tr.xx = []; tr.ii = []; tr.ll = {}; op = [];
                %extract user-defined inputs
                if (nargin >= +2), op = varargin{1}; end
                %user-defined inputs
                if (~isstruct(op))
                    op.nobj = +32;
                    op.long = .75;
                else
                %bound population count
                    if (isfield(op,'nobj'))
                        if (op.nobj <= +0 )
                            error('maketree:invalidInputs','Invalid options.');
                        end
                    else
                        op.nobj = +32;
                    end
                %bound "long" tolerance
                    if (isfield(op,'long'))
                        if (op.long < +.0 || op.long > +1.)
                            error('maketree:invalidInputs','Invalid options.');
                        end
                    else
                        op.long = .75; 
                    end
                end

            %dimensions of rectangles
            nd = size(rp,2) / +2 ;
            ni = size(rp,1) ;
            %alloc. workspace
            xx = zeros(ni*1,2*nd);
            ii = zeros(ni*1,2);
            ll = cell (ni*1,1);
            ss = zeros(ni*1,1);
            %min & max coord. masks
            lv = false(size(rp,2),1);
            rv = false(size(rp,2),1);
            lv((1:nd)+nd*+0) = true ;
            rv((1:nd)+nd*+1) = true ;

            %rectangle centres
            rc =(rp(:,lv)+rp(:,rv)) * +.5;
            %rectangle lengths
            rd = rp(:,rv)-rp(:,lv);

            %indexing for root node
            ii(1,1) = +0;
            ii(1,2) = +0;
            %root contains all rectangles
            ll{1}    = (+1:ni)' ;
            xx(1,lv) = min(rp(:,lv),[],1);
            xx(1,rv) = max(rp(:,rv),[],1);

            %divide nodes until all constraints satisfied
            ss(+1) = +1; ns = +1; nn = +1;   
            while (ns ~= +0)
                %pop node from stack
                ni = ss(ns); 
                ns = ns - 1;
                %divide node if too populous
                if (length(ll{ni}) > op.nobj)
                %child indexing
                    n1 = nn+1; 
                    n2 = nn+2;
                %set of rectangles in parent
                    li = ll{ni} ;
                %split plane on longest axis
                    ndim = xx(ni,rv)-xx(ni,lv);
                   [mx,ax] = max(ndim);

                %push rectangles to children
                    il = rd(li,ax) > op.long * mx;
                    lp = li( il);           %  "long" rectangles
                    ls = li(~il);           % "short" rectangles

                %select split position: take the mean of the set of
                %(non-"long") rectangle centres along axis AX

                    sp = sum(rc(ls,ax))/length(ls);
                %partition based on centres
                    i2 = rc(ls,ax)>sp ;
                    l1 = ls(~i2);           %  "left" rectangles
                    l2 = ls( i2);           % "right" rectangles

                %deal with empty partitions
                    if (isempty(l1) || isempty(l2))
                        continue;
                    end

                %finalise node position
                    xx(n1,lv) = min(rp(l1,lv),[],1);
                    xx(n1,rv) = max(rp(l1,rv),[],1);
                    xx(n2,lv) = min(rp(l2,lv),[],1);
                    xx(n2,rv) = max(rp(l2,rv),[],1);
                %parent--child indexing
                    ii(n1,1) = ni; 
                    ii(n2,1) = ni;
                    ii(ni,2) = n1;
                %finalise list indexing
                    ll{ni}   = lp;
                    ll{n1}   = l1;
                    ll{n2}   = l2;
                %push child nodes onto stack
                    ss(ns+1) = n1;
                    ss(ns+2) = n2;
                    
                    nn = nn+2; ns = ns+2;
                end       
            end
            
            xx           = xx(1:nn,:);
            ii           = ii(1:nn,:);
            ll(nn+1:end) = [] ;
            
            tr.xx = xx; 
            tr.ii = ii; 
            tr.ll = ll;
            
        end
        
        function [tm,im] = scantree(tr,pi,fn)
            
            tm.ii = []; tm.ll = {};
            im.ii = []; im.ll = {};
            
            tm.ii = zeros(size(tr.ii,1),1);
            tm.ll = cell (size(tr.ii,1),1); 
            im.ii = zeros(size(pi,1),1);
            im.ll = cell (size(pi,1),1);
            ss    = zeros(size(tr.ii,1),1);
            sl    = cell (size(tr.ii,1),1);
            sl{1} = (1:size(pi,1))';
    
            tf = ~cellfun('isempty',tr.ll);

            %partition query items amongst tree nodes: from root
            ss(1) = +1; ns = +1; no = +1;
            while (ns ~= +0)
                %pop node from stack
                ni = ss(ns); ns = ns - 1;

                if (tf(ni))
                    %push output: non-empty node NI contains items LL
                    tm.ii(no) = ni; 
                    tm.ll{no} = sl{ns+1};
                    no = no + 1 ;
                end

                if (tr.ii(ni,+2)~=+0)
                    %partition amongst child nodes
                    c1 = tr.ii(ni,+2)+0;
                    c2 = tr.ii(ni,+2)+1;
                    %user-defined partitions of LL
                   [j1,j2] = feval(fn,pi(sl{ns+1},:), ...
                                   tr.xx(c1,:), ...
                                   tr.xx(c2,:)) ;
                    %lists of items per child node
                    l1 = sl{ns+1}(j1); 
                    l2 = sl{ns+1}(j2);

                    if (~isempty(l1))
                        %push nonempty node onto stack
                        ns = ns + 1; 
                        ss(ns) = c1; 
                        sl{ns} = l1; 
                    end
                    if (~isempty(l2)) 
                        %push nonempty node onto stack
                        ns = ns + 1; 
                        ss(ns) = c2; 
                        sl{ns} = l2; 
                    end
                end

            end
                
            tm.ii(no:end) = [];
            tm.ll(no:end) = [];

            %compute inverse map only if desired
            if (nargout == +1), return; end

            %accumulate pair'd tree-item matches
            ic = cell(no-1,+1); jc = tm.ll;
            for ip = +1 : no-1
                ni = tm.ii(ip);
                ic{ip} = ni(ones(length(jc{ip}),1));
            end
            ii = vertcat(ic{:}); ni = size(ii,1);
            jj = vertcat(jc{:}); nj = size(jj,1);
            %invert ordering via sort
            [jj,jx] = sort (jj) ; 
            ii = ii(jx);
            jx = find(diff(jj)~=+0);
            im.ii = [jj(jx);jj(nj)];
            jx = [+0;jx;ni];
            %distribute single item-tree matches
            for ip = +1 : size(im.ii,1)
                im.ll{ip} = ii(jx(ip+0)+1:jx(ip+1));
            end
                
        end
                
        function [ip,it] = testpts(pp,tt,pi,pk,tk)
        %TESTPTS compute the point-tria matches within a tile.

            mp = length(pk); 
            mt = length(tk);

            switch (size(tt,2))
                case 3
                %pts in 2-simplexes
                pk = pk';
                pk = pk(ones(mt,1),:); pk = pk(:);
                tk = tk(:,ones(1,mp)); tk = tk(:);

                in = intria2(pp,tt(tk,:),pi(pk,:));
                ip = pk(in);
                it = tk(in);

                case 4
                %pts in 3-simplexes
                pk = pk';
                pk = pk(ones(mt,1),:); pk = pk(:);
                tk = tk(:,ones(1,mp)); tk = tk(:);

                in = intria3(pp,tt(tk,:),pi(pk,:));
                ip = pk(in);
                it = tk(in);
            end

        end

        function [in] = intria2(pp,tt,pi)
        %INTRIA2 returns TRUE for points enclosed by 2-simplexes.

            vi = pp(tt(:,1),:)-pi;
            vj = pp(tt(:,2),:)-pi;
            vk = pp(tt(:,3),:)-pi;
            %compute sub-volume about PI
            aa = zeros(size(tt,1),3) ;
            aa(:,1) = vi(:,1).*vj(:,2) - ...
                      vj(:,1).*vi(:,2) ;
            aa(:,2) = vj(:,1).*vk(:,2) - ...
                      vk(:,1).*vj(:,2) ;
            aa(:,3) = vk(:,1).*vi(:,2) - ...
                      vi(:,1).*vk(:,2) ;
            %PI is internal if same sign
        
            in = all(aa>=0-tol,2) | all(aa<=0+tol,2);
      

        end

        function [in] = intria3(pp,tt,pi)
        %INTRIA3 returns TRUE for points enclosed by 3-simplexes.

            v1 = pi-pp(tt(:,1),:);
            v2 = pi-pp(tt(:,2),:);
            v3 = pi-pp(tt(:,3),:);
            v4 = pi-pp(tt(:,4),:);
            %compute sub-volume about PI
            aa = zeros(size(tt,1),4) ;
            aa(:,1) = ...
            +v1(:,1).*(v2(:,2).*v3(:,3) - ...
                      v2(:,3).*v3(:,2) ) ...
            -v1(:,2).*(v2(:,1).*v3(:,3) - ...
                      v2(:,3).*v3(:,1) ) ...
            +v1(:,3).*(v2(:,1).*v3(:,2) - ...
                      v2(:,2).*v3(:,1) ) ;
            aa(:,2) = ...
            +v1(:,1).*(v4(:,2).*v2(:,3) - ...
                      v4(:,3).*v2(:,2) ) ...
            -v1(:,2).*(v4(:,1).*v2(:,3) - ...
                      v4(:,3).*v2(:,1) ) ...
            +v1(:,3).*(v4(:,1).*v2(:,2) - ...
                      v4(:,2).*v2(:,1) ) ;
            aa(:,3) = ...
            +v2(:,1).*(v4(:,2).*v3(:,3) - ...
                      v4(:,3).*v3(:,2) ) ...
            -v2(:,2).*(v4(:,1).*v3(:,3) - ...
                      v4(:,3).*v3(:,1) ) ...
            +v2(:,3).*(v4(:,1).*v3(:,2) - ...
                      v4(:,2).*v3(:,1) ) ;
            aa(:,4) = ...
            +v3(:,1).*(v4(:,2).*v1(:,3) - ...
                      v4(:,3).*v1(:,2) ) ...
            -v3(:,2).*(v4(:,1).*v1(:,3) - ...
                      v4(:,3).*v1(:,1) ) ...
            +v3(:,3).*(v4(:,1).*v1(:,2) - ...
                      v4(:,2).*v1(:,1) ) ;
      
            in = all(aa>=0-tol,2) | all(aa<=0+tol,2);


        end
        
        xq1(:,1,:) = xq'; %reshape querypoint-structure to compute baryzentrics
% 		idx = elidx(idx);
        bary = obj.barycentric_coordinates(xq1, idx);
        
        end 
        
        function bary = barycentric_coordinates(obj,x,idx)
            % BARYCENTRIC_COORDINATES returns the barycentric coordinates
            % of points with coordinates x located in elements of index
            % idx. x is suposed to be a matrixarray with the same structure
            % as mesh.co. idx is an indexvector containing the element's
            % numbers x is in.
            %
            obj.jacobiandata();
            bary = obj.DinvT(:,:,idx)'*(x-obj.co(:,1,obj.el(idx,1)));
            bary = [1-sum(bary,1); bary];
            bary = double(permute(bary,[3,1,2]));
        end
     
    end
end








