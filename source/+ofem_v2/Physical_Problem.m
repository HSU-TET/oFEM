classdef Physical_Problem < handle
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        element;
        geometry;
        
        has_stiffness;
        paraS;
        S;
        SPart = {};
        
        has_damping;
        D;
        DPart = {};
        
        has_mass;
        paraM;
        M;
        MPart = {};
        
        F = {};
        b;
        
        M_robin;
        
        bc;
        ic;
        
        DOFs;
        
        u;
        
        matlib;
    end
    
    methods
        
        
        function obj = Physical_Problem(element, geometry, has_s, has_m, has_d)
            obj.element = element;
            obj.geometry = geometry;
            obj.has_stiffness = has_s;
            obj.has_damping = has_d;
            obj.has_mass = has_m;
            obj.u = zeros(obj.geometry.Nco,1);
            
        end
        
        function setParaS(obj, namePara)
            obj.paraS = namePara;
        end
        
        function setParaM(obj,namePara)
            obj.paraM = namePara;
        end
        
        function setBoundaryCondition(obj, BC)
            n= size(obj.geometry.bd, 2);
            obj.geometry.bd{3,BC.index} = BC;
        end
        
        function assemble(obj)
            %Nc = size(obj.geometry.co  ,3);
            Np = size(obj.geometry.parts, 2);
            
            %obj.DOFs = 1:Nc;
            Nc = obj.DOFs.Nd;
            
            obj.S=sparse(Nc,Nc);
            obj.D=sparse(Nc,Nc);
            obj.M=sparse(Nc,Nc);
            obj.b=sparse(Nc,1);
            
            obj.u = sparse(Nc,1);
            
            obj.M_robin=sparse(Nc,Nc);
            
            has_f = size(obj.geometry.parts,1)>=4;
            
            % start assembling
            tic;
            
            for i=1:Np
                pIdx = obj.geometry.parts{3,i};
                if obj.has_stiffness
                    obj.SPart{i} = obj.element.assembleStiffness(obj,pIdx,obj.geometry.parts{2,i}.(obj.paraS));
                    obj.S = obj.S + obj.SPart{i};
                end
                if obj.has_mass
                    obj.MPart{i} = obj.element.assembleMass(obj,pIdx,obj.geometry.parts{2,i}.(obj.paraM));
                    obj.M = obj.M+obj.MPart{i};
                end
                if has_f
                    if ~isempty(obj.geometry.parts{4,i})
                        obj.F{i} = obj.geometry.parts{4,i}.force(obj,pIdx);
                        obj.b = obj.b+obj.F{i};
                    end
                end                
            end
            
            if size(obj.geometry.bd,1) == 3
                nB= size(obj.geometry.bd,2);
                for i=1:nB
                    if ~isempty(obj.geometry.bd{3,i})
                        obj.u = obj.u + obj.geometry.bd{3,i}.loadVector(obj);
                    end
                end
                obj.b = obj.b - (obj.S+obj.M+obj.D)*obj.u;
            end
        end
        
        function attachDOFHandler(obj,DOFHandler)
            obj.DOFs = DOFHandler;
        end
        
        function du = gradCell(obj,u)
            % Computes the gradient inside each cell. No patch recovery
            DinvT = obj.geometry.DinvT;
            uElem = u(obj.geometry.el(:,:));
            uElem = ofem_v2.tools.matrixarray(reshape(uElem',size(uElem,2),1,[]));
            phi = obj.element.gradShapeFunctions;
            du = (DinvT*phi')*uElem;
        end
        
        
        
        function solve(obj)
            %DOFs = obj.DOFs.getDOFs;
            DOFs = obj.DOFs;
            obj.u(DOFs) = obj.S(DOFs,DOFs)\obj.b(DOFs);
            
        end
        
        function parabolic_solve(obj, init, timesteps, deltaTime)
            
            %           Passt das so ?berhaupt?!
            
            obj.u(obj.dof)= init.value;
            dt= deltaTime;
            n = timesteps;
            for i = 0:n
                b2 = obj.b+obj.M/dt*obj.u;
                obj.u(obj.dof) = (obj.S(obj.dof,obj.dof)+obj.M(obj.dof, obj.dof)/dt+obj.M_robin(obj.dof, obj.dof))\b2(obj.dof); % change something
                
                obj.geometry.export_UCD([pwd, '/export'],['export',num2str(i)], {'T', obj.u, ''});
            end
        end
        
    end
end

