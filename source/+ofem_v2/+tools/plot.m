function [h] = plot(obj,u,varargin)
%This function plots the solution u of the created Geometry class object.
%Varargin can be used to label the figures legend and z-axis in 2D.


co = double(squeeze(obj.co)');

        name = 'u';
        if nargin==3
            if ~ischar(varargin{1})
                warning('ofem:plot:InvalidArgument',...
                        'Optional argument must be a valid string!');
            else
                name=varargin{1};
            end
        end

        switch obj.dim
            case 2
                h=trimesh(obj.el          , ...
                          co(:,1)              , ...
                          co(:,2)              , ...
                          full(u)                    , ...
                          'FaceColor', 'interp', ...
                          'EdgeColor', 'none'  );
                xlabel('x-axis');
                ylabel('y-axis');
                zlabel(name);
                legend(name);
                colorbar;

            case 3
               warning('ofem:plot:3D solution plot is not supported. Use export_UCD instead.');
        end


end

