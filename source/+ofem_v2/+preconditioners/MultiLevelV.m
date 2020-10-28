classdef MultiLevelV
    %V-Cycle multilevel preconditioner
    %Start at level 1, go to top, then descend to 1 again. Not the other
    %way around
    
    properties
        dofs = {};
        levels = 0;
    end
    
    methods 
        function obj = MultiLevelV(fe,dofs)
            obj.dofs{1} = dofs.e2DOF(:,1);
            obj.levels = fe.degree;
            obj.dofs{1} = setdiff(obj.dofs{1},dofs.fixedDOFs);
            
            if obj.levels < 1
                disp('why tho?')
            else
                for i = 1:obj.levels
                    obj.dofs{i+1} = unique(dofs.e2DOF(:,i+1));
                    if i > 1
                        obj.dofs{i+1} = [obj.dofs{i+1};unique(dofs.f2DOF(:,(i-2)*3+1:(i-1)*3))];
                        for j = i+1:obj.levels
                            obj.dofs{i+1} = [obj.dofs{i+1};dofs.f2DOF(:,(levels-1)*3+1+j*2:(levels-1)*3+1+j*2)];
                        end
                    end
                    obj.dofs{i} = setdiff(obj.dofs{i},dofs.fixedDOFs);
                end
            end
        end
        
        function x = apply(obj,A,b,fdofs,varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if nargin > 4
                S = varargin{1};
                M = varargin{2};
            end
            a = zeros(size(A,1),1);
            if length(a)~=length(b)
                a(fdofs) = b;
            else
                a = b;
            end
            if exist('S')
                x = obj.multilevel(A,a,S,M);
            else
                x = obj.multilevel(A,a);
            end
            if length(x) ~= length(b)
                x = x(fdofs);
            end
        end
    
        function e = multilevel(obj,A,b,varargin)
            if nargin > 3
                S = varargin{1};
                M = varargin{2};
%                 L = ichol(S(obj.dofs{1},obj.dofs{1})+M(obj.dofs{1},obj.dofs{1}));
            end
            e = zeros(size(b));
            p = obj.levels;
            % upcycle
            for m = 1:obj.levels+1
                if m == 1
                    e(obj.dofs{m}) = A(obj.dofs{m},obj.dofs{m})\b(obj.dofs{m});
%                     [e(obj.dofs{m}),~] = bicgstab(A(obj.dofs{m},obj.dofs{m}),b(obj.dofs{m}),1e-4,100,L,L',e(obj.dofs{m}));
                else
%                     e(obj.dofs{m}) = obj.SSOR(A(obj.dofs{m},obj.dofs{m}),b(obj.dofs{m})); %presmooth
                    e(obj.dofs{m}) = obj.SSOR(A(obj.dofs{m},obj.dofs{m}),b(obj.dofs{m})); % postsmooth
                    %[e(obj.dofs{m}),~] = bicgstab(A(obj.dofs{m},obj.dofs{m}),b(obj.dofs{m}),1e-2);
                end
                for n = m+1:p+1
                    b(obj.dofs{n}) = b(obj.dofs{n}) - A(obj.dofs{n},obj.dofs{m})*e(obj.dofs{m});
                end
            end
            % downcycle
            for m = p:-1:1
                for n = m+1:p
                    b(obj.dofs{m}) = b(obj.dofs{m}) - A(obj.dofs{m},obj.dofs{n})*e(obj.dofs{n});
                end
                if m == 1
                    e(obj.dofs{m}) = A(obj.dofs{m},obj.dofs{m})\b(obj.dofs{m});
%                     [e(obj.dofs{m}),~] = bicgstab(A(obj.dofs{m},obj.dofs{m}),b(obj.dofs{m}),1e-4,100,L,L',e(obj.dofs{m}));
                else
%                     e(obj.dofs{m}) = obj.SSOR(A(obj.dofs{m},obj.dofs{m}),b(obj.dofs{m})); % postsmooth
                    e(obj.dofs{m}) = obj.SSOR(A(obj.dofs{m},obj.dofs{m}),b(obj.dofs{m})); % postsmooth
                    %[e(obj.dofs{m}),~] = bicgstab(A(obj.dofs{m},obj.dofs{m}),b(obj.dofs{m}),1e-2);
                end
            end
%             if i == 0
%                 x(obj.dofs{i+1}) = A(obj.dofs{i+1},obj.dofs{i+1})\b(obj.dofs{i+1});
%             else
%                 x(obj.dofs{i+1}) = SSOR(A(obj.dofs{i+1},obj.dofs{i+1}),b(obj.dofs{i+1}),1e-2,10);
%                 b(obj.dofs{i}) = b(obj.dofs{i}) - A(obj.dofs{i},obj.dofs{i+1})*x(obj.dofs{i+1});
%                 x = obj.multilevel(A,b,i-1);
%                 b(obj.dofs{i+1}) = b(obj.dofs{i+1}) - A(obj.dofs{i+1},obj.dofs{i})*x(obj.dofs{i});
%                 x(obj.dofs{i+1}) = SSOR(A(obj.dofs{i+1},obj.dofs{i+1}),b(obj.dofs{i+1}),1e-2,10);
%             end   
        end
        
        function x = decompType(obj,A,b)
            L = ichol(sqrt(A.*conj(A)));
            z = L\b;
            x = L'\z;
        end
        
        function x = SSOR(obj,A,b)
            L = tril(A,-1);
            U = triu(A,1);
            D = diag(diag(A));
            w = 1.2;
            x = zeros(size(b));
            %M = (1/(2-w))*(1/w*D+L)*(1/w*D)\(1/w*D+L).';
%             M = (1/w*D+L).';
%             M = (1/w*D)\M;
%             M = (1/w*D+L)*M;
%             M = (1/(2-w))*M;
            %M = D+w*L;
            %N = D+w*U;
            K = sqrt(2-w)*diag((sqrt(diag(D/w))).^(-1))*(eye(size(D))-L*diag(diag(D/w).^(-1)));
            
            for i = 1:5
                %z = (eye(size(D))-D\L)*(b-A*x);
                %x = x + (eye(size(D))-L/D)*z;
                x = x + K'*(K*(b-A*x));
            end
        end
        
        function x = gaussSeidel(obj,A,b)
            L = tril(A);
            U = triu(A,1);
            x = zeros(size(b));
            
            for i = 1:5
                x = -L\U*x + L*b;
            end
        end
    end
end
























