classdef edgefx
    %EDGEFX Build a mesh size function gridded interpolant from a velocity
    %       model
    
    %-----------------------------------------------------------
    %   Keith Roberts   : 2019 --
    %   Email           : krober@usp.br
    %   Last updated    : 10/20/2019
    %-----------------------------------------------------------
    %
    
    properties(Access=private)
        wl  % size of elements (m) per speed of sound in sea water.
        hvp % mesh size as a function of seismic p-wavespeed.
        used % edge function keywords
        feat % geodata class instance
        min_el % minimum mesh size (m)
        max_el % maximum mesh size (m)
        g % mesh size gradation rate (decimal percent)
    end
    
    properties(Access=public)
        F % gridded interpolant mesh size function
    end
    
    
    methods(Access=public)
        function obj = edgefx(varargin)
            %
            p = inputParser;
            
            defval = 0; % placeholder value if arg is not passed.
            
            % add name/value pairs
            addOptional(p,'wl',defval);
            addOptional(p,'min_el',defval);
            addOptional(p,'max_el',defval);
            addOptional(p,'g',0.35);
            addOptional(p,'geodata',defval);
            
            % parse the inputs
            parse(p,varargin{:});
            % store the inputs as a struct
            inp = p.Results;
            % get the fieldnames of the edge functions
            inp = orderfields(inp,{'geodata','wl','min_el','max_el','g'});
            flds = fieldnames(inp);
            for i = 1 : numel(flds)
                type = flds{i};
                switch type
                    case('wl')
                        obj.wl = inp.(flds{i});
                        if obj.wl~=0
                            obj.wl = inp.(flds{i});
                        end
                        assert(obj.wl > 0); 
                    case('min_el')
                        obj.min_el = inp.(flds{i});
                        if obj.min_el~=0
                            obj.min_el = inp.(flds{i});
                        end
                        assert(obj.min_el > 0);
                    case('max_el')
                        obj.max_el = inp.(flds{i}); 
                        if obj.max_el~=0
                            obj.max_el = inp.(flds{i}); 
                        end
                        assert(obj.max_el > obj.min_el);
                    case('geodata')
                        if isa(inp.(flds{i}),'geodata')
                            obj.feat = inp.(flds{i});
                        else
                            error('GEODATA OBJECT REQUIRED');
                        end
                    case('g')
                        obj.g = inp.(flds{i});
                        if obj.g~=0.35 
                            obj.g = inp.(flds{i});
                        end
                        assert(obj.g < 1);
                        assert(obj.g > 0);
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
        
    end % end non-static public methods
    
    methods(Access=private)
        %% Wavelength edgefx.
        % speed of sound in sea water gets wl element size (in m).
        function obj = wlfx(obj)
            vp_sw = 1484; % m/s speed of sound in sea water
            
            [yg,zg] = obj.feat.CreateStructGrid;
            Fvp = GetFvp(obj.feat);
            tmp = Fvp(yg,zg);
            obj.hvp=(obj.wl*tmp)./vp_sw;
        end
        
        
        function obj = finalize(obj)
            ny = obj.feat.GetNy; 
            nz = obj.feat.GetNz; 
            gsp = obj.feat.GetGridspace; 
            
            disp('Enforcing size bounds...'); 
            
            hh = obj.hvp ; 
            hh(hh<obj.min_el)=obj.min_el;
             hh(hh>obj.max_el)=obj.max_el;

            disp('Relaxing the mesh size gradient...');
            hfun = zeros(ny*nz,1);
            nn = 0;
            for ipos = 1 : ny
                for jpos = 1 : nz
                    nn = nn + 1;
                    hfun(nn,1) = hh(ipos,jpos);
                end
            end
            dy = repmat(gsp,[1,nz]); % for gradient function
            dz = gsp;        % for gradient function
            [hfun,flag] = obj.limgradStruct(nz,dy,dz,hfun,...
                obj.g,sqrt(length(hfun)));
            if flag == 1
                disp('Gradient relaxing converged!');
            else
                error(['FATAL: Gradient relaxing did not converge, '
                    'please check your edge functions']);
            end
            % reshape it back
            nn = 0;
            for ipos = 1 : ny
                for jpos = 1 : nz
                    nn = nn+1;
                    hh(ipos,jpos) = hfun(nn);
                end
            end
            clearvars hfun fdfdx
            
            % 
            hh_m=hh;
            
            [xg,yg]=obj.feat.CreateStructGrid;
            obj.F = griddedInterpolant(xg,yg,hh_m,'linear','nearest');
            clearvars xg yg

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
    
