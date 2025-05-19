classdef Coercivity < ofem_v2.Volume_force
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Properties
    
    %% Methods
    methods
        
        % function obj= Coercivity(value, part, mesh)
        %     if isscalar(part)
        %         partIndex = part;
        %     else
        %         [~, number_colums] = size(mesh.parts);
        %         for i= 1:number_colums
        %             if isequal(mesh.parts{1,i}, part)
        %                 partIndex = i;
        %                 break;
        %             end
        %         end
        %     end
        %     if isa(value,'ofem_v2.tools.matrixarray')
        %         obj.value = value(:,:,mesh.parts{3, partIndex});
        %     else
        %         obj.value = value;
        %     end
        %     obj.part = partIndex;
        %     mesh.parts{4, partIndex} = obj;
        % end

        function f = force(obj, phys, pIdx)
            f = phys.element.derivativeForce(phys,obj.value,pIdx);
        end
                
    end
end


















