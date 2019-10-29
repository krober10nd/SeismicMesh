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
                                   'min_el','max_el','g','cr','dt'});
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
                        if obj.f~=10
                            obj.f = inp.(flds{i});
                            assert(obj.f > 0);
                        end
                    case('slp')
                        obj.slp = inp.(flds{i});
                        if obj.slp~=0
                            obj.slp = inp.(flds{i});
                            assert(obj.f > 0);
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
            [yg,zg] = obj.feat.CreateStructGrid;
            Fvp = GetFvp(obj.feat);
            vp = Fvp(yg,zg);
            tmp = stdfilt(vp); % uses a default filter of 3x3
            tmp = (1-tmp./max(tmp(:)));
            obj.hslp = tmp*obj.slp + obj.min_el;
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
                hh = zeros([nx,ny,nz,1])+obj.min_el;;
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
                 hfun = reshape(hh_m',[numel(hh_m),1]);
                 dy = repmat(gsp,[1,nz]); % for gradient function
                 dz = gsp;        % for gradient function
                 [hfun,flag] = obj.limgradStruct(nz,dy,dz,hfun,...
                     obj.g,sqrt(length(hfun)));
                 assert(flag ==1);
                 disp('Gradient relaxing converged!');
                 % reshape it back
                 nn = 0;
                 for ipos = 1 : ny
                     for jpos = 1 : nz
                         nn = nn+1;
                         hh_m(ipos,jpos) = hfun(nn);
                     end
                 end
                 clearvars hfun fdfdx
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
   function [ffun,flag] = limgradStruct(ny,xeglen,yeglen,ffun,fdfdx,imax)
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
            
            if length(fdfdx)==1
                fdfdx = ffun*0 + fdfdx ;
            end
            %----------------------------- ASET=ITER if node is "active"
            aset = zeros(size(ffun,1),1) ;
            
            %----------------------------- exhaustive 'til all satisfied
            ftol = min(ffun) * sqrt(eps) ;
            
            % ------------------------ calculate hypotenuse eglen length
            dgeglen = sqrt(xeglen.^2+repmat(yeglen,1,ny).^2);
            eglen   = 0.5*(xeglen+repmat(yeglen,1,ny));
            
            rm = zeros(9,1);
            for iter = +1 : imax
                
                %------------------------- find "active" nodes this pass
                aidx = find(aset == iter - 1) ;
                
                if (isempty(aidx)), break; end
                %------------------------- reorder => better convergence
                [aval,idxx] = sort(ffun(aidx)) ;
                
                aidx = aidx(idxx);
                
                %------------------------- speed up a little by preallocating occasionally
                npos = zeros(9,1); elens = zeros(9,1);
                
                %------------------------- visit adj. edges and set DFDX
                for i = 1 : length(aidx)
                    % ----- map doubly index to singly indexed
                    inod = aidx(i);
                    ipos = 1 + floor((inod-1)/ny);
                    jpos = inod - (ipos - 1)*ny;
                    
                    % ------ use 8 edge stencil for higher order gradients
                    nn=1;
                    npos(nn) =  inod;                        nn=nn+1;
                    npos(nn) =  ipos*ny    + jpos;           nn=nn+1;%--- nnod of right adj
                    npos(nn) = (ipos-2)*ny + jpos;           nn=nn+1;%--- nnod of left adj
                    npos(nn) = (ipos-1)*ny + min(jpos+1,ny); nn=nn+1;%--- nnod of above adj
                    npos(nn) = (ipos-1)*ny + max(jpos-1,1);  nn=nn+1;%--- nnod of below adj
                    npos(nn) = (ipos*ny)   +  max(jpos-1,1); nn=nn+1;%--- nnod of right bot diag adj
                    npos(nn) =  (ipos*ny)  + min(jpos+1,ny); nn=nn+1;%--- nnod of right top diag adj
                    npos(nn) = (ipos-2)*ny +  min(jpos+1,ny);nn=nn+1;%--- nnod of left top diag adj
                    npos(nn) = (ipos-2)*ny + max(jpos-1,1);          %--- nnod of left bot diag adj
                    
                    
                    % ------ populate elens
                    nn = 1;
                    elens(nn) = eglen(jpos);  nn=nn+1;
                    elens(nn) = xeglen(jpos); nn=nn+1;
                    elens(nn) = xeglen(jpos); nn=nn+1;
                    elens(nn) = yeglen; nn=nn+1;
                    elens(nn) = yeglen; nn=nn+1;
                    elens(nn) = dgeglen(jpos); nn=nn+1;
                    elens(nn) = dgeglen(jpos); nn=nn+1;
                    elens(nn) = dgeglen(jpos); nn=nn+1;
                    elens(nn) = dgeglen(jpos);
                    
                    %----- handle boundary vertex adjs.
                    rm = npos <= 0 | npos > size(ffun,1);
                    npos(rm) = [];
                    elens(rm)= [];
                    
                    for ne = 2 : length(npos)
                        nod1 = npos(1);
                        nod2 = npos(ne);
                        elen = elens(ne);
                        
                        %----------------- calc. limits about min.-value
                        if (ffun(nod1) > ffun(nod2))
                            
                            fun1 = ffun(nod2) ...
                                + elen * fdfdx(nod2) ;
                            
                            if (ffun(nod1) > fun1+ftol)
                                ffun(nod1) = fun1;
                                aset(nod1) = iter;
                            end
                            
                        else
                            
                            fun2 = ffun(nod1) ...
                                + elen * fdfdx(nod2) ;
                            
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
    
