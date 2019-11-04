classdef edgefx
    %EDGEFX Build a mesh size function gridded interpolant from a velocity
    %       model
    
    %-----------------------------------------------------------
    %   Keith Roberts   : 2019 --
    %   Email           : krober@usp.br
    %   Last updated    : 10/27/2019
    %-----------------------------------------------------------
    %
    
    properties(Access=private)
        wl  % number of nodes given a wave with frequency f 
        f   % maximum frequency of source injection
        hvp % mesh size as a function of seismic p-wavespeed.
        slp % number of nodes per gradient 
        hslp % mesh size as a function of seismic p-wavespeed gradient
        used % edge function keywords
        feat % geodata class instance
        min_el % minimum mesh size (m)
        max_el % maximum mesh size (m)
        g % mesh size gradation rate (decimal percent)
        cr % courant number 
        dt % desired timestep 
    end
    
    properties(Access=public)
        F % gridded interpolant mesh size function
    end
    
    
    methods(Access=public)
        function obj = edgefx(varargin)
            %
            p = inputParser;
            
            % add name/value pairs
            addOptional(p,'wl',0);
            addOptional(p,'f',0); 
            addOptional(p,'slp',0); 
            addOptional(p,'min_el',0);
            addOptional(p,'max_el',inf);
            addOptional(p,'g',0);
            addOptional(p,'dt',0); 
            addOptional(p,'cr',0); 
            addOptional(p,'geodata',0);
            
            % parse the inputs
            parse(p,varargin{:});
            % store the inputs as a struct
            inp = p.Results;
            % get the fieldnames of the edge functions
            inp = orderfields(inp,{'geodata','wl','f',...
                                   'slp',...
                                   'min_el','max_el',...
                                   'g','cr','dt'});
            flds = fieldnames(inp);
            for i = 1 : numel(flds)
                type = flds{i};
                switch type
                    case('wl')
                        obj.wl = inp.(flds{i});
                        if obj.wl~=0
                            obj.wl = inp.(flds{i});
                            assert(obj.wl > 0);
                        end
                    case('f')
                        obj.f = inp.(flds{i});
                        if obj.f~=0
                            obj.f = inp.(flds{i});
                            assert(obj.f > 0);
                        end
                    case('slp')
                        obj.slp = inp.(flds{i});
                        if obj.slp~=0
                            obj.slp = inp.(flds{i});
                            assert(obj.slp > 0);
                        end
                    case('min_el')
                        obj.min_el = inp.(flds{i});
                        if obj.min_el~=0
                            obj.min_el = inp.(flds{i});
                            assert(obj.min_el > 0);
                        end
                    case('max_el')
                        obj.max_el = inp.(flds{i});
                        if ~isinf(obj.max_el)
                            obj.max_el = inp.(flds{i});
                            assert(obj.max_el > obj.min_el);
                        end
                    case('geodata')
                        if isa(inp.(flds{i}),'geodata')
                            obj.feat = inp.(flds{i});
                        else
                            error('GEODATA OBJECT REQUIRED');
                        end
                    case('g')
                        obj.g = inp.(flds{i});
                        if obj.g~=0
                            assert(obj.g < 1);
                        end
                    case('cr')
                        obj.cr = inp.(flds{i}); 
                        if obj.cr~=0
                          assert(obj.cr < 1); 
                        end
                    case('dt')
                        obj.dt = inp.(flds{i});
                        if obj.dt~=0
                          assert(obj.dt > 0);
                          assert(obj.cr > 0); % require cr to be set if dt is set.
                        end
                end
            end
            
            
            % now turn on the edge functions
            for i = 1 : numel(flds)
                type = flds{i};
                switch type
                    case('wl')
                        obj.wl  = inp.(flds{i});
                        if obj.wl(1)~=0
                            disp('Building wavelength function...');
                            obj = wlfx(obj);
                            obj.used{end+1} = 'wl';
                        end
                    case('slp')
                        obj.slp  = inp.(flds{i});
                        if obj.slp(1)~=0
                            disp('Building slope function...');
                            obj = slpfx(obj);
                            obj.used{end+1} = 'slp';
                        end
                end
            end
            
            obj = finalize(obj);
            
        end
        
        
        function [axH]=plot(obj)
            [yg,zg]=obj.feat.CreateStructGrid ;
            skip=5 ; % save memory and time by skipping
            figure;
            axH=pcolor(yg(1:skip:end,1:skip:end),...
                zg(1:skip:end,1:skip:end),...
                obj.F.Values(1:skip:end,1:skip:end)) ;
            shading interp;
            set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
            xlabel('Y-position (m)');
            ylabel('Z-position/depth (m)');
            cb=colorbar; ylabel(cb,'Desired mesh size (m)') ;
            set(gca,'FontSize',16) ;
            axis equal; axis tight; 
        end
        
        function min_el = GetMinEl(obj); min_el=obj.min_el; end
        
        function max_el = GetMaxEl(obj); max_el=obj.max_el; end 
        
        function feat = GetFeat(obj); feat=obj.feat; end 

    end % end non-static public methods
    
    methods(Access=private)
        %% Wavelength edgefx.
        function obj = wlfx(obj)
            
            if obj.feat.GetDim == 2
                [yg,zg] = obj.feat.CreateStructGrid;
                Fvp = GetFvp(obj.feat);
                vp = Fvp(yg,zg);
                obj.hvp=(vp.*obj.f)./obj.wl;
            else
                [xg,yg,zg] = obj.feat.CreateStructGrid3D;
                Fvp = GetFvp(obj.feat);
                vp = Fvp(xg,yg,zg);
                obj.hvp=(vp.*obj.f)./obj.wl;
            end
        end
        %% slope edge function
        function obj = slpfx(obj)
            if obj.feat.GetDim == 2
                [yg,zg] = obj.feat.CreateStructGrid;
                Fvp = GetFvp(obj.feat);
                vp = Fvp(yg,zg);
                tmp = stdfilter(vp,[3,3],1,'replicate'); % use a larger stencil here for smoother variance
                tmp = real(tmp);
                tmp = (1-tmp./500); %  demoninator is reference gradient for which min_el is mapped
                obj.hslp = max(tmp*obj.slp + obj.min_el,10); % ensure result doesn't go negative.
            else
                [xg,yg,zg] = obj.feat.CreateStructGrid3D;
                Fvp = GetFvp(obj.feat);
                vp = Fvp(xg,yg,zg);
                for k = 1 : size(zg,1)
                    tmp = stdfilter(squeeze(vp(k,:,:)),[3,3],1,'replicate'); % use a larger stencil here for smoother variance
                    tmp = real(tmp);
                    tmp = (1-tmp./500); %  demoninator is reference gradient for which min_el is mapped
                    obj.hslp(k,:,:) = max(tmp*obj.slp + obj.min_el,10); % ensure result doesn't go negative.
                end
            end
        end

        
        function obj = finalize(obj)
            
            nx = obj.feat.GetNx; 
            ny = obj.feat.GetNy; 
            nz = obj.feat.GetNz; 
            gsp = obj.feat.GetGridspace; 
            dim = obj.feat.GetDim; 
            
            % package the edge functions into a known order.
            counter = 0;
            if dim==2
                hh = zeros([ny,nz,1])+obj.min_el;
            else
                hh = zeros([nx,ny,nz,1])+obj.min_el;
            end
            for i = 1 : numel(obj.used)
                type = obj.used{i};
                switch type
                    case('wl')
                        counter = counter + 1;
                        if dim==2
                            hh(:,:,counter) = obj.hvp;
                        else
                            hh(:,:,:,counter) = obj.hvp;
                        end
                        obj.hvp = single(obj.hvp);
                    case('slp')
                        counter = counter + 1;
                        if dim==2
                            hh(:,:,counter) = obj.hslp;
                        else
                            hh(:,:,:,counter) = obj.hslp;
                        end
                        obj.hslp = single(obj.hslp);
                    otherwise
                        error('FATAL:  Could not finalize edge function');
                end
            end
            
            if dim==2
                [hh_m] = min(hh,[],3);
            else
                [hh_m] = min(hh,[],4);
            end
            clearvars hh
            
             disp('Enforcing size bounds...');
             hh_m(hh_m<obj.min_el)=obj.min_el;
             hh_m(hh_m>obj.max_el)=obj.max_el;
             
             if obj.g > 0
                 disp('Relaxing the mesh size gradient...');
                 hfun = reshape(hh_m,[numel(hh_m),1]); % reshape column wise
                 hfun = FastHJ( int32([ny nz 1]), gsp, obj.g, int32(sqrt(length(hfun))), hfun);
                 disp('Gradient relaxing converged!');
                 hh_m = reshape(hfun,[ny,nz]);
                 clearvars hfun
             end
             
             % enforce the CFL if present
             if obj.dt > 0
                 disp(['Enforcing timestep of ',num2str(obj.dt),' seconds.']);
                 if dim==2
                     [yg,zg] = obj.feat.CreateStructGrid;
                     Fvp = GetFvp(obj.feat);
                     vp = Fvp(yg,zg);
                 else
                     [xg,yg,zg] = obj.feat.CreateStructGrid3D;
                     Fvp = GetFvp(obj.feat);
                     vp = Fvp(xg,yg,zg);
                 end
                 %  vp*dt/dx < cr (generally < 1 for stability).
                 cfl = (obj.dt*vp)./hh_m; % this is your cfl
                 dxn = vp.*obj.dt/obj.cr;      % assume simulation time step of dt sec and cfl of dcfl;
                 hh_m( cfl > obj.cr) = dxn( cfl > obj.cr);   %--in planar metres
                 clear cfl dxn u hh_d;
             end
            
            if dim==2
                [yg,zg]=obj.feat.CreateStructGrid;
                obj.F = griddedInterpolant(yg,zg,hh_m,'linear','nearest');
                clearvars yg zg
            else
                [xg,yg,zg]=obj.feat.CreateStructGrid3D;
                obj.F = griddedInterpolant(xg,yg,zg,hh_m,'linear','nearest');
                clearvars xg yg zg
            end

        end
                
    end % end private non-static methods
    
    methods(Static,Access=private)
        % Suppose you have a 3D matrix A,
        % and you want to get the singleton index for A(i1, i2, i3)
        % which is given by  sub2ind(size(A), i1, i2, i3).
        % This is the equivalent expression:
        % i1 + (i2-1)*size(A,1) + (i3-1)*size(A,1)*size(A,2)
        
        function [ffun,flag] = limgradStruct(ny,elen,ffun,fdfdx,imax)
            %LIMGRAD impose "gradient-limits" on a function defined over
            %an undirected graph.
            %         Modified for a structred graph with eight node stencil
            %         Last updated: 24/06/2017
            %         Email       : krober10@nd.edu
            %         Keith Roberts, 2017.
            % ---------------------
            %             Modified to have spatially variable fdfdx
            %             Last updated: 27/04/2019
            %             Keith Roberts, 2019
            
            %----------------------------- ASET=ITER if node is "active"
            aset = zeros(size(ffun,1),1) ;
            
            %----------------------------- exhaustive until all are satisfied
            ftol = min(ffun) * sqrt(eps) ;
            
            rm = zeros(5,1);
            
            % -----
            elenXfdfdx = elen*fdfdx; 
            
            for iter = 1 : imax
                
                %------------------------- find "active" nodes this pass
                aidx = find(aset == iter - 1) ;
                
                if (isempty(aidx)), break; end
                %------------------------- reorder => better convergence
                [~,idxx] = sort(ffun(aidx)) ;
                
                aidx = aidx(idxx);
                
                %------------------------- speed up a little by preallocating occasionally
                npos = zeros(5,1);
                
                %------------------------- 
                for i = 1 : length(aidx)
                                        
                    % ----- map doubly index to singly indexed
                    inod = aidx(i);
                    ipos = 1 + floor((inod-1)/ny);
                    jpos = inod - (ipos - 1)*ny;

                    
                    % ------ gather indices use 4 edge stencil
                    npos(1) =  inod;                       
                    npos(2) =  ipos*ny    + jpos;           %--- nnod of right adj
                    npos(3) = (ipos-2)*ny + jpos;           %--- nnod of left adj
                    npos(4) = (ipos-1)*ny + min(jpos+1,ny); %--- nnod of above adj
                    npos(5) = (ipos-1)*ny + max(jpos-1,1);  %--- nnod of below adj
                    
                    %----- handle boundary vertex adjs.
                    rm = npos <= 0 | npos > size(ffun,1);
                    npos(rm) = [];
                                            
                    nod1 = npos(1);
                    ffunNod1 = ffun(nod1); 
                    
                    for ne = 2 : length(npos)
                        nod2 = npos(ne);                        
                        
                        %----------------- calc. limits about min.-value
                        if ( ffunNod1 > ffun(nod2))
                            
                            fun1 = ffun(nod2) ...
                                + elenXfdfdx;
                            
                            if (ffunNod1 > fun1+ftol)
                                ffunNod1 = fun1;
                                aset(nod1) = iter;
                            end
                            
                        else
                            
                            fun2 = ffunNod1 ...
                                + elenXfdfdx ;
                            
                            if (ffun(nod2) > fun2+ftol)
                                ffun(nod2) = fun2;
                                aset(nod2) = iter;
                            end
                            
                        end
                    end
                end
                rm = rm*0;
                flag = (iter < imax) ;
                
            end
        end % end limgradstruct 
    end %% end static methods
end %% end class
    
