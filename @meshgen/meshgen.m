classdef meshgen < msh
    %   MESHGEN: Mesh generation class
    %   Handles input parameters to create a meshgen class object that can be
    %   used to build a msh class.
    
    properties(Access=private)
        feat % geodata object
        egfx % edgefx object
        anno % Approx. Nearest Neighbor search object.
        annData % datat contained with KD-tree in anno
        plot_on % to plot, or to not plot, that is the question
    end
    
    methods(Access=private)
        function obj = meshgen(varargin)
            p = inputParser;
            % unpack options and set default ones, catch errors.
            
            defval = 0; % placeholder value if arg is not passed.
            
            addOptional(p,'plot_on',defval);
            addOptional(p,'edgefx',defval);
            addOptional(p,'geodata',defval);
            
            % parse the inputs
            parse(p,varargin{:});
            
            inp = p.Results;
            
            inp = orderfields(inp,{ });
            % get the fieldnames of the edge functions
            fields = fieldnames(inp);
            
            for i = 1 : numel(fields)
                type = fields{i};
                switch type
                    case('plot_on')
                        obj.plot_on= inp.(fields{i});
                    case('edgefx')
                        if isa(inp.(fields{i}),'edgefx')
                            obj.feat = inp.(fields{i});
                        else
                            error('GEODATA OBJECT REQUIRED');
                        end
                    case('geodata')
                        if isa(inp.(fields{i}),'geodata')
                            obj.egfx = inp.(fields{i});
                        else
                            error('EDGEFX OBJECT REQUIRED');
                        end
                end
            end
        end
    end % end prviate non-state members
    
    methods(Access=public)
        function  obj = build(obj)
            %DISTMESH2D 2-D Mesh Generator using Distance Functions.
            %%
            % FORM INITIAL POINTS HERE
            
            
            N = size(p,1); % Number of points N
            disp(['Number of initial points after rejection is ',num2str(N)]);
            %% Iterate
            pold = inf;                                                    % For first iteration
            if obj.plot_on >= 1
                clf,view(2),axis equal;
            end
            toc
            fprintf(1,' ------------------------------------------------------->\n') ;
            disp('Begin iterating...');
            while 1
                tic
                if ~mod(it,obj.nscreen)
                    disp(['Iteration =' num2str(it)]) ;
                end
                
                % 3. Retriangulation by the Delaunay algorithm
                if max(sqrt(sum((p(1:size(pold,1),:)-pold).^2,2))/h0) > ttol         % Any large movement?
                    p = fixmesh(p);                                        % Ensure only unique points.
                    N = size(p,1); pold = p;                               % Save current positions
                    [t,p] = obj.delaunay_elim(p,obj.fd,geps,0);                % Delaunay with elimination
                    
                    N = size(p,1);
                    
                    % 4. Describe each bar by a unique pair of nodes.
                    bars = [t(:,[1,2]); t(:,[1,3]); t(:,[2,3])];           % Interior bars duplicated
                    bars = unique(sort(bars,2),'rows');                    % Bars as node pairs
                    
                    % 5. Graphical output of the current mesh
                    if obj.plot_on >= 1 && (mod(it,obj.nscreen)==0 || it == 1)
                        cla,triplot(p(:,1),p(:,2),t)
                        title(['Iteration = ',num2str(it)]);
                    end
                end
                
                % Getting element quality and check goodness
                
                % Termination quality, mesh quality reached is copacetic.
                if ~mod(it,imp)
                    qual_diff = mq_l3sig - obj.qual(max(1,it-imp),2);
                    if abs(qual_diff) < obj.qual_tol
                        % Do the final elimination of small connectivity
                        [t,p] = delaunay_elim(p,obj.fd,geps,1);
                        disp('Quality of mesh is good enough, exit')
                        close all;
                        break;
                    end
                end
                
                % Saving a temp mesh
                if ~mod(it,obj.nscreen) && delIT == 0
                    disp(['Number of nodes is ' num2str(length(p))])
                    disp(['Mean mesh quality is ' num2str(mq_m)])
                    disp(['Min mesh quality is ' num2str(mq_l)])
                    disp(['3rd sigma lower mesh quality is ' num2str(mq_l3sig)])
                    tempp = p; tempt = t;
                    save('Temp_grid.mat','it','tempp','tempt');
                    clearvars tempp tempt
                end
                
                % 6. Move mesh points based on bar lengths L and forces F
                barvec = pt(bars(:,1),:)- pt(bars(:,2),:);                 % List of bar vectors
                L
                ideal_bars = (p(bars(:,1),:) + p(bars(:,2),:))/2;          % Used to determine what bars are in bbox
                hbars = 0*ideal_bars(:,1);
                
                
                L0 = hbars*Fscale*median(L)/median(hbars);                  % L0 = Desired lengths using ratio of medians scale factor
                LN = L./L0;                                                 % LN = Normalized bar lengths
                
                % Mesh improvements (deleting and addition)
                if ~mod(it,imp)
                    nn = []; pst = [];
                    if qual_diff < imp*0.01 && qual_diff > 0
                        % Remove elements with small connectivity
                        nn = get_small_connectivity(p,t);
                        disp(['Deleting ' num2str(length(nn)) ' due to small connectivity'])
                        
                        % Remove points that are too close (< LN = 0.5)
                        if any(LN < 0.5)
                            % Do not delete pfix too close.
                            nn1 = setdiff(reshape(bars(LN < 0.5,:),[],1),[(1:nfix)']);
                            disp(['Deleting ' num2str(length(nn1)) ' points too close together'])
                            nn = unique([nn; nn1]);
                        end
                        
                        % Split long edges however many times to
                        % better lead to LN of 1
                        if any(LN > 2)
                            nsplit = floor(LN);
                            nsplit(nsplit < 1) = 1;
                            adding = 0;
                            % Probably we can just split once?
                            for jj = 2:2
                                il = find(nsplit >= jj);
                                xadd = zeros(length(il),jj-1);
                                yadd = zeros(length(il),jj-1);
                                for jjj = 1 : length(il)
                                    deltax = (p(bars(il(jjj),2),1)- p(bars(il(jjj),1),1))/jj;
                                    deltay = (p(bars(il(jjj),2),2)- p(bars(il(jjj),1),2))/jj;
                                    xadd(jjj,:) = p(bars(il(jjj),1),1) + (1:jj-1)*deltax;
                                    yadd(jjj,:) = p(bars(il(jjj),1),2) + (1:jj-1)*deltay;
                                end
                                pst = [pst; xadd(:) yadd(:)];
                                adding = numel(xadd) + adding;
                            end
                            disp(['Adding ',num2str(adding) ,' points.'])
                        end
                    end
                    if ~isempty(nn) || ~isempty(pst)
                        % Doing the actual subtracting and add
                        p(nn,:)= [];
                        p = [p; pst];
                        pold = inf;
                        it = it + 1;
                        continue;
                    end
                end
                
                F    = (1-LN.^4).*exp(-LN.^4)./LN;                         % Bessens-Heckbert edge force
                Fvec = F*[1,1].*barvec;
                
                Ftot = full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
                Ftot(1:nfix,:) = 0;                                        % Force = 0 at fixed points
                pt = pt + deltat*Ftot;                                     % Update node positions
                
                [p(:,1),p(:,2)] = m_xy2ll(pt(:,1),pt(:,2));
                
                %7. Bring outside points back to the boundary
                d = feval(obj.fd,p,obj,[],1); ix = d > 0;                  % Find points outside (d>0)
                ix(1:nfix)=0;
                if sum(ix) > 0
                    dgradx = (feval(obj.fd,[p(ix,1)+deps,p(ix,2)],obj,[])...%,1)...
                        -d(ix))/deps; % Numerical
                    dgrady = (feval(obj.fd,[p(ix,1),p(ix,2)+deps],obj,[])...%,1)...
                        -d(ix))/deps; % gradient
                    dgrad2 = dgradx.^+2 + dgrady.^+2;
                    p(ix,:) = p(ix,:)-[d(ix).*dgradx./dgrad2,...
                        d(ix).*dgrady./dgrad2];
                end
                
                % 8. Termination criterion: Exceed itmax
                it = it + 1 ;
                
                if ( it > obj.itmax )
                    % Do the final deletion of small connectivity
                    [t,p] = delaunay_elim(p,obj.fd,geps,1);
                    disp('Max iterations reached, exiting')
                    close all;
                    break ;
                end
                toc
            end
            %%
            disp('Finished iterating...');
            fprintf(1,' ------------------------------------------------------->\n') ;
            
            %% Doing the final cleaning and fixing to the mesh...
            % Clean up the mesh if specified
            if ~strcmp(obj.cleanup,'none')
                % Put the mesh class into the grd part of meshgen and clean
                obj.grd.p = p; obj.grd.t = t;
                [obj.grd,qout] = clean(obj.grd,obj.cleanup,...
                    'nscreen',obj.nscreen,'djc',obj.dj_cutoff,...
                    'pfix',obj.pfix);
                obj.grd.pfix = obj.pfix ;
                obj.grd.egfix= obj.egfix ;
                obj.qual(end+1,:) = qout;
            else
                % Fix mesh on the projected space
                [p(:,1),p(:,2)] = m_ll2xy(p(:,1),p(:,2));
                [p,t] = fixmesh(p,t);
                [p(:,1),p(:,2)] = m_xy2ll(p(:,1),p(:,2));
                % Put the mesh class into the grd part of meshgen
                obj.grd.p = p; obj.grd.t = t;
                obj.grd.pfix = obj.pfix ;
                obj.grd.egfix= obj.egfix ;
            end
            
            
        end
    end
    
    
    methods(static,Access=private)
        
        function [t,p] = delaunay_elim(p,fd,geps,final)
            % Removing mean to reduce the magnitude of the points to
            % help the convex calc
            if exist('pt1','var'); clear pt1; end
            [pt1(:,1),pt1(:,2)] = m_ll2xy(p(:,1),p(:,2));
            if isempty(obj.egfix)
                p_s  = pt1 - repmat(mean(pt1),[N,1]);
                TR   = delaunayTriangulation(p_s);
            else
                TR   = delaunayTriangulation(pt1(:,1),pt1(:,2),obj.egfix);
                pt1  = TR.Points;
            end
            for kk = 1:final+1
                if kk > 1
                    % Perform the following below upon exit from the mesh
                    % generation algorithm
                    nn = get_small_connectivity(pt1,t);
                    nn1 = heal_fixed_edges(pt1,t,obj.egfix) ;
                    nn = unique([nn; nn1]) ;
                    TR.Points(nn,:) = [];
                    pt1(nn,:) = [];
                end
                t = TR.ConnectivityList;
                pmid = squeeze(mean(reshape(pt1(t,:),[],3,2),2));      % Compute centroids
                [pmid(:,1),pmid(:,2)] = m_xy2ll(pmid(:,1),pmid(:,2));  % Change back to lat lon
                t    = t(feval(fd,pmid,obj,[]) < -geps,:);             % Keep interior triangles
                if kk == 1
                    % Deleting very straight triangles
                    tq_n = gettrimeshquan( pt1, t);
                    bad_ele = any(tq_n.vang < 1*pi/180 | ...
                        tq_n.vang > 179*pi/180,2);
                    t(bad_ele,:) = [];
                end
            end
            if length(pt1) ~= length(p)
                clear p
                [p(:,1),p(:,2)] = m_xy2ll(pt1(:,1),pt1(:,2));
            end
        end
        
        function nn = get_small_connectivity(p,t)
            % Get node connectivity (look for 4)
            [~, enum] = VertToEle(t);
            % Make sure they are not boundary nodes
            bdbars = extdom_edges2(t, p);
            bdnodes = unique(bdbars(:));
            I = find(enum <= 4);
            nn = setdiff(I',[(1:nfix)';bdnodes]);
            return;
        end
        
    end % end static-private methods
    
    
end % end meshgen class