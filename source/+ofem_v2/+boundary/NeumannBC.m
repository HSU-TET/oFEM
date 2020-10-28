classdef Neumann < handle
    
    %% Proteries
     properties
         value; 
         boundary; 
         index; 
         meas; 
         normalVector; 
		 bd;
         
         bNeumann;
     end
     
     %% Methods
     methods
         
         function obj=Neumann(value, boundaryName, mesh)
            obj.value = value;   
             
            n = size(mesh.bd, 2);
            for i = 1:n
                if isequal(mesh.bd{1,i},boundaryName)
					obj.bd = mesh.bd{2,i};
                    obj.boundary = unique(mesh.bd{2,i});
                    obj.index= n; 
                    break;
                end
            end  
            
            switch mesh.type
                case 'tri'
                    obj.neumann1Tri(mesh);
             
                case 'tet'
                    obj.neumann1Tet(mesh);
                    
                otherwise
                    error('ofem:mesh:Unspecified',...
                          'Unspecified error found');   
            end
         end
         
         function obj = loadVector(obj, physicalProblem)
        %pressure returns the pressure-originated part of the load vector.
        %
        % b=pressure(meas,phi,w,l,g,faces,co) computes the pressure
        % originated vector in terms of a quadrature rule. w and l are the
        % weights and quadrature points of the rule, respectively. phi
        % contains the values of the shape functions evaluated at the
        % quadrature points and g is a functions handle returning the
        % pressure at arbitrary points.
        
%         l = physicalProblem.element.gausianLength; 
%         w = physicalProblem.element.gausianWeight; 
%         phi = physicalProblem.element.shape_function; 
%         co = physicalProblem.geometry.co; 
%         g = obj.value; 
%         meas=?
%         faces?!

	        w = [1; 1]/2;
            l = 0.5+[-1; 1]/(2*sqrt(3));
            l = [1-l, l];
			
			shape = l';
               
            Nl     = size(l  ,2); % number of barycentric coordinates
            Nc     = size(physicalProblem.geometry.co ,3);
            Nq     = size(w  ,1);
            Ns     = size(shape,1);
            Nf     = size(obj.bd,1);

            % faceco gives global quadrature points => eval there
            faceco = reshape(physicalProblem.geometry.co(:,:,obj.bd(:,1:Nl)'),[],Nl,Nf);

            F      = ofem_v2.tools.matrixarray(zeros(1,Ns,Nf));

            for q=1:Nq
                X = faceco*(l');
                F = F + obj.value*(w(q)*shape(:,q)');
            end

            F = F*obj.meas;
            obj.boundary = obj.bd';
            obj.bNeumann = sparse(obj.boundary(:),1,F(:),Nc,1);
         end
         
         function neumann1Tri(obj, mesh)
             % calculates and assigns meas and normalVector to Properties
             %
             edges = mesh.co(:,:,mesh.bd{2, obj.index}(:,2))- mesh.co(:,:,mesh.bd{2, obj.index}(:,1));
             obj.meas = sqrt(dot(edges,edges,1));
             obj.normalVector = -edges.rot;
         end
         
         function neumann1Tet(obj, mesh)
             % calculates and assigns meas and normalVector to Properties
             %
            e1         = mesh.co(:,:,mesh.bd{2, obj.index}(:,2))-mesh.co(:,:,mesh.bd{2, obj.index}(:,1));
            e2         = obj.co(:,:,obj.boundary(:,3))-obj.co(:,:,obj.boundary(:,1));
            obj.normalVector = cross(e1,e2,1);
            obj.meas    = sqrt(dot(obj.normalVector,obj.normalVector,1))/2;
         end
         
        function [meas,faces,normals,opp]=neumann(obj, physicalProblem)
        %NEUMANN returns the normals associated to the Neumann boundary.
        %
        % data=NEUMANN(idx), where idx is an index vector indexing the bd
        % cell, returns a cell array data the same length as idx. In each
        % cell of data we have the following data per row:
        % [ length, normal, IDs ],
        % where length denotes the length of the respective edge on the
        % reference element, normal is the normals vector of the global
        % edge and IDs contains the indices of the edges' nodes.
        % Each cell of data corresponds to a requested sideset.
        %
            ss      = obj.bd(2,idx); %boundaries
            Nss     = length(idx);
            meas    = cell(Nss);
            opp     = cell(Nss);
            normals = cell(Nss);
            faces   = cell(Nss);
            npere   = size(obj.el,2);

            switch obj.type
                case 'edge'
                    %% edges

                case 'tri'
                    switch obj.filetype
                                       
                        case 'msh'
                            for i=1:Nss
                                switch npere
                                    case 3
                                        %% first order mesh
                                        % 1 2 3
                                        faces{i} = ss{i};

                                    case 6
                                        %% second order mesh
                                        % 1 2 3 e12 e23 e31
                                        faces{i} = [obj.el(ss{1,i}{2,1},[1 2 4]);... % E1
                                                    obj.el(ss{1,i}{2,2},[2 3 5]);... % E2
                                                    obj.el(ss{1,i}{2,3},[3 1 6])];   % E3

                                end

                                edges      = obj.co(:,:,faces{i}(:,2))-obj.co(:,:,faces{i}(:,1));
                                meas{i}    = sqrt(dot(edges,edges,1));
                                normals{i} = -edges.rot;


                            end
                    end
                            
                case 'quad'
                    %% quadrilateral
                    error('ofem:mesh:InvalidMesh',...
                          'Quadrilateral meshes not supported so far!');
                case 'tet'
                    %% tets
                    switch obj.filetype
                        case 'inp'
                            % sideset coding for 3D tetrahedral structures by
                            % the Abaqus inp file is as follows:
                            % S1: 1, 2, 3
                            % S2: 1, 2, 4
                            % S3: 2, 3, 4
                            % S4: 1, 3, 4
                            for i=1:Nss
                                switch npere
                                    case 4
                                        %% first order mesh
                                        % 1 2 3 4
                                        faces{i} = ss;

                                    case 10
                                        %% second order mesh
                                        % 1 2 3 4 e12 e13 e14 e23 e24 e34
                                        faces{i} = [obj.el(ss{1,i}{2,1},[1 2 3 5  8 6]);... % S1
                                                    obj.el(ss{1,i}{2,2},[1 2 4 5  9 7]);... % S2
                                                    obj.el(ss{1,i}{2,3},[2 3 4 8 10 9]);... % S3
                                                    obj.el(ss{1,i}{2,4},[1 3 4 6 10 7])];   % S4

                                end

                                opp{i}   = [ss(:,4); ... % S1
                                            ss(:,3); ... % S2
                                            ss(:,1); ... % S3
                                            ss(:,2)];    % S4

                                e1         = obj.co(:,:,faces{i}(:,2))-obj.co(:,:,faces{i}(:,1));
                                e2         = obj.co(:,:,faces{i}(:,3))-obj.co(:,:,faces{i}(:,1));
                                normals{i} = cross(e1,e2,1);
                                meas{i}    = sqrt(dot(normals{i},normals{i},1))/2;


                                % correct direction
                                tmp = ofem_v2.tools.matrixarray(ones(1,1,size(normals{i},3)));
                                tmp(1,1,dot(normals{i},obj.co(:,:,opp{i}),1)>=0) = -1;

                                normals{i} = normals{i}./(2*meas{i}.*tmp);
                            end
                        case 'msh'

                            for i=1:Nss
                                switch npere
                                    case 4
                                        %% first order mesh
                                        % 1 2 3 4
                                        faces{i} = ss{i};

                                    case 10
                                        %% second order mesh
                                        % 1 2 3 4 e12 e13 e14 e23 e24 e34
                                        faces{i} = [obj.el(ss{1,i}{2,1},[1 2 3 5  8 6]);... % S1
                                                    obj.el(ss{1,i}{2,2},[1 2 4 5  9 7]);... % S2
                                                    obj.el(ss{1,i}{2,3},[2 3 4 8 10 9]);... % S3
                                                    obj.el(ss{1,i}{2,4},[1 3 4 6 10 7])];   % S4

                                end

%                                 opp{i}   = [obj.el(ss{1,i}{2,1},4); ... % S1
%                                             obj.el(ss{1,i}{2,2},3); ... % S2
%                                             obj.el(ss{1,i}{2,3},1); ... % S3
%                                             obj.el(ss{1,i}{2,4},2)];    % S4

                                e1         = obj.co(:,:,faces{i}(:,2))-obj.co(:,:,faces{i}(:,1));
                                e2         = obj.co(:,:,faces{i}(:,3))-obj.co(:,:,faces{i}(:,1));
                                normals{i} = cross(e1,e2,1);
                                meas{i}    = sqrt(dot(normals{i},normals{i},1))/2;


                                % correct direction
%                                 tmp = ofem.matrixarray(ones(1,1,size(normals{i},3)));
%                                 tmp(1,1,dot(normals{i},obj.co(:,:,opp{i}),1)>=0) = -1;
% 
%                                 normals{i} = normals{i}./(2*meas{i}.*tmp);
                            end
                    end
 
                case 'hex'
                    % _________________hexahedron__________________________
                    error('ofem:mesh:InvalidMesh',...
                          'Hexahedral meshes not supported so far!');
                otherwise
                    error('ofem:mesh:Unspecified',...
                          'Unspecified error found');
            end
        end
     
     end
    
    
end