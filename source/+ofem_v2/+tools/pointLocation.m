function [idx, tr, bary] = pointLocation(obj, xq, opt)
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
%This function is a modification of https://github.com/dengwirda/find-tria


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
            tt  = [tt;obj.el(obj.parts{3,opt.parts(i)},:)]; %#ok<AGROW>
            elidx = [elidx,obj.parts{3,opt.parts(i)}]; %#ok<AGROW>
       end
    else
       tt = obj.el;
    end
    
    pp  = double(reshape(permute(obj.co,[3,1,2]),[],size(obj.co,1)));
    tol = 1e-12;
    tp  = [];  %#ok<NASGU>
    tj  = [];  %#ok<NASGU>
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

    ic     = cell(size(tm.ii,1),1);
    jc     = cell(size(tm.ii,1),1);
    
    for ip = 1 : size(tm.ii,1) %#ok<*FXUP>
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

    xq1(:,1,:) = xq'; %reshape querypoint-structure to compute baryzentrics
    % 		idx = elidx(idx);
    bary = obj.barycentric_coordinates(xq1, idx);

    %% functions
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

    function [tr] = maketree(rp,varargin) 

        tr.xx = []; tr.ii = []; tr.ll = {}; op = [];

        %extract user-defined inputs
        if (nargin >= +2)
            op = varargin{1}; 
        end

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
                li = ll{ni};

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
                               tr.xx(c2,:));
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
        if (nargout == +1)
            return; 
        end

        %accumulate pair'd tree-item matches
        ic = cell(no-1,+1); 
        jc = tm.ll;

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
        aa = zeros(size(tt,1),3);

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

end