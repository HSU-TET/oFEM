function [meshOut] = shrinkMesh(meshIn,parts,varargin)
%SHRINKMESH Summary of this function goes here
%   Detailed explanation goes here

	meshOut = ofem_v2.Geometry;
	idx = [];
	j = 1;
	for i = parts
		idx = [idx;meshIn.parts{3,i}(:)];
		meshOut.parts{1,j} = meshIn.parts{1,i};
		meshOut.parts{2,j} = meshIn.parts{2,i};
		meshOut.parts{3,j} = meshIn.parts{3,i};
		j = j+1;
	end

	
end

