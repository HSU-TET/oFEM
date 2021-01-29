classdef Neumann < handle
	
	%% Proteries
	properties
		value;
		boundary;
		index;
		meas;
		normalVector;
		bd;
		dim;
		
		b;
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
					edges = mesh.co(:,:,mesh.bd{2, obj.index}(:,2))- mesh.co(:,:,mesh.bd{2, obj.index}(:,1));
					obj.meas = sqrt(dot(edges,edges,1));
					obj.normalVector = -edges.rot;
					obj.dim = 2;
					
				case 'tet'
					e1 = mesh.co(:,:,mesh.bd{2, obj.index}(:,2))-mesh.co(:,:,mesh.bd{2, obj.index}(:,1));
					e2 = mesh.co(:,:,obj.boundary(:,3))-mesh.co(:,:,obj.boundary(:,1));
					obj.normalVector = cross(e1,e2,1);
					obj.meas = sqrt(dot(obj.normalVector,obj.normalVector,1))/2;
					obj.dim = 3;
					
				otherwise
					error('ofem:mesh:Unspecified',...
						'Unspecified error found');
			end
		end
		
		function obj = loadVector(obj, physicalProblem)			
			[w,l] = ofem_v2.tools.gaussSimplex(obj.dim-1, 1);
			
			
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
	end
	
end