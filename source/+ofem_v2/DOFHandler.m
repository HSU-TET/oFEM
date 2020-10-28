classdef DOFHandler < handle
	%DOFHANDLER What does the DOF handler handle? The DOF handler handles
	%DOFs!
	%   Detailed explanation goes here
	
	properties
		DOFs;
		nDOFs;
        n2DOF;
		eDOFs;
		e2DOF;
		fDOFs;
		f2DOF;
		iDOFs;
		i2DOF;
        el2DOF;
        Nd;
		boundaryDOFs;
		elements = {};
		geom;
		freeDOFs;
        fixedDOFs = [];
	end
	
	methods
		function obj = DOFHandler(geom)
			%DOFHANDLER Construct an instance of this class
			%   Detailed explanation goes here
			obj.geom = geom;
		end
		
		function attach(obj,elem)
			%METHOD1 Summary of this method goes here
			%   Detailed explanation goes here
			obj.elements{size(obj.elements,2)+1} = elem;
		end
		
		function generateDOFs(obj)
			obj.nDOFs = 1:obj.elements{1}.nodeDOFs*obj.geom.Nco;
			obj.eDOFs = 1:obj.elements{1}.edgeDOFs/6*obj.geom.Ned;
            obj.eDOFs = reshape(obj.eDOFs,[],obj.elements{1}.edgeDOFs/6)';
            obj.eDOFs = obj.eDOFs(:)';
			obj.fDOFs = 1:obj.elements{1}.faceDOFs/4*obj.geom.Nfa;
			obj.iDOFs = 1:obj.elements{1}.interiorDOFs*obj.geom.Nint;
			
			if ~isempty(obj.eDOFs)&&~isempty(obj.nDOFs)
				obj.eDOFs = obj.eDOFs+max(obj.nDOFs);
			end
			if ~isempty(obj.fDOFs)&&~isempty(obj.eDOFs)
				obj.fDOFs = obj.fDOFs+max(obj.eDOFs);
			end
			if ~isempty(obj.iDOFs)&&~isempty(obj.fDOFs)
				obj.iDOFs = obj.iDOFs+max(obj.fDOFs);
            elseif ~isempty(obj.iDOFs) && ~isempty(obj.eDOFs)
                obj.iDOFs = obj.iDOFs+max(obj.eDOFs);
            end
            if ~isempty(obj.nDOFs)
                obj.n2DOF = reshape(obj.nDOFs,obj.elements{1}.nodeDOFs,[])';
                obj.el2DOF = [obj.el2DOF,...
                    reshape(obj.n2DOF(obj.geom.el',:)',[],obj.geom.Nint)'];
            end
            if ~isempty(obj.eDOFs)
                obj.e2DOF = reshape(obj.eDOFs,obj.elements{1}.edgeDOFs/6,[])';
                obj.el2DOF = [obj.el2DOF,...
                    reshape(obj.e2DOF(obj.geom.el2ed(:),:),obj.geom.Nint,obj.elements{1}.edgeDOFs)];
            end
            if ~isempty(obj.fDOFs)
                obj.f2DOF = reshape(obj.fDOFs,obj.elements{1}.faceDOFs/4,[])';
                obj.el2DOF = [obj.el2DOF,...
                    reshape(obj.f2DOF(obj.geom.el2fa',:)',[],obj.geom.Nint)'];
            end
            if ~isempty(obj.iDOFs)
                obj.i2DOF = reshape(obj.iDOFs,obj.elements{1}.interiorDOFs,[])';
                obj.el2DOF = [obj.el2DOF,obj.i2DOF];
            end
			obj.DOFs = [obj.nDOFs,obj.eDOFs,obj.fDOFs,obj.iDOFs];
            obj.Nd = max(obj.DOFs);
			obj.freeDOFs = obj.DOFs;
		end
		
		function DOFs = getDOFs(obj)
			DOFs = obj.freeDOFs;
        end
        
        function DOFs = AMSDOFs(obj,phys)
            Nco = phys.geometry.Nco;
            Ned = phys.geometry.Ned;
            Nfa = phys.geometry.Nfa;
            deg = phys.element.degree;
            nDofs = 1:Nco;
            eDofs = [];
            fDofs = [];
            iDofs = [];
            for i = 1:size(phys.geometry.bd,2)
                nDofs = setdiff(nDofs,phys.geometry.bd{3,i}.nodes);
                if deg > 0
                    eDofs = 1:Ned;
                    eDofs = setdiff(eDofs,phys.geometry.bd{3,i}.edges);
                    eDofs = eDofs'+((0:deg-1)*Ned);
                    eDofs = sort(eDofs(:))+Nco;
                end
                if deg > 1 && phys.element.faceDOFs > 0
                    fDofs = 1:Nfa;
                    fDofs = setdiff(fDofs,phys.geometry.bd{3,i}.faces);
                    fDofs = obj.f2DOF(fDofs,:);
                    fDofs = sort(fDofs(:))-Ned+Nco;
                end
                if deg > 2 && phys.element.interiorDOFs > 0
                    iDofs = unique(obj.iDOFs(:));
                    iDofs = iDofs + Nco-Ned;
                end
            end
            DOFs = [nDofs,eDofs',fDofs',iDofs'];
        end
		
		function reduceDOFs(obj,BD)
            obj.boundaryDOFs{BD.index} = [];
			if ~isempty(BD.nodes)&&~isempty(obj.nDOFs)
                n = obj.n2DOF(BD.nodes,:);
                obj.boundaryDOFs{BD.index} = [obj.boundaryDOFs{BD.index};n(:)];
                obj.freeDOFs = setdiff(obj.freeDOFs,n(:));
            end
            if ~isempty(BD.edges)&&~isempty(obj.eDOFs)
                e = obj.e2DOF(BD.edges,:);
                obj.boundaryDOFs{BD.index} = [obj.boundaryDOFs{BD.index};e(:)];
                obj.freeDOFs = setdiff(obj.freeDOFs,e(:));
            end
            if ~isempty(BD.faces)&&~isempty(obj.fDOFs)
                f = obj.f2DOF(BD.faces,:);
                obj.boundaryDOFs{BD.index} = [obj.boundaryDOFs{BD.index};f(:)];
                obj.freeDOFs = setdiff(obj.freeDOFs,f(:));
            end
            if ~isempty(BD.interiors)&&~isempty(obj.iDOFs)
                i = obj.i2DOF(BD.interiors,:);
                obj.boundaryDOFs{BD.index} = [obj.boundaryDOFs{BD.index};i(:)];
                obj.freeDOFs = setdiff(obj.freeDOFs,i(:));
            end
            obj.boundaryDOFs{BD.index} = unique(obj.boundaryDOFs{BD.index});
            obj.fixedDOFs = [obj.fixedDOFs;obj.boundaryDOFs{BD.index}];
		end
		
	end
end




















