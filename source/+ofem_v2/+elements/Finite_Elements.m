classdef (Abstract) Finite_Elements < handle
    % Abstract class of different finite element types
    %   Detailed explanation goes here
    %% Properties
    
    properties (Abstract, Access = public)
        
        % Shape function of the element
        shape_function;
    end
    %% Methods
    
    methods (Abstract)
        computeBasis(shape_function)
		assembleMass(shape_function)
		assembleStiffness(shape_function)
		assembleDamping(shape_function)
    end
    
    methods
        
       function [co,el]=create_midpoints(obj)
            Nel=numel(obj.el(:,1));
            
            co = obj.co;

            switch obj.type
                case 'edge'
                    % edge
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
                    % triangle
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
                    % quadrilateral
                    error('ofem:mesh:NotImplemented',...
                          'Quadrilateral meshes not supported so far!');
                    
                case 'tet'
                    % tetrahedron
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
                    % hexahedron
                    error('ofem:mesh:NotImplemented',...
                          'Hexahedral meshes not supported so far!');
                    
                otherwise
                    error('ofem:mesh:Unspecified',...
                          'Unspecified error found');
            end
         end
      
    end
end


