classdef geodata
    %GEODATA Read in a velocity model from a SEGY file
    
    %-----------------------------------------------------------
    %   Keith Roberts   : 2019 --
    %   Email           : krober@usp.br
    %   Last updated    : 10/28/2019
    %-----------------------------------------------------------
    %
    properties(Access=private)
        fname % path to segy seismic data.
        x0y0  % top left corner coordinate
        Fvp   % gridded interpolant of seismic wavespeed in km/s
        ny    % number of grid points in y-direction
        nx    % number of grid points in the x-direction
        nz    % number of grid points in z-direction
        gridspace    % grid space (in m)
        dim   % dimension of problem 
    end
    
    properties(Access=public)
        bbox % domain ranges in x y (and z)
        expand % amount of expand the domain in the -x, +x, and -z directions
    end
        
    methods(Access=public)
        % default class constructor
        % GEODATA construct the default class.
        % loads in a velocity model from a file.
        function obj = geodata(varargin)
            %
            p = inputParser;
            
            defval=0;
            
            % add name/value pairs
            addOptional(p,'velocity_model',defval);
            addOptional(p,'expand',0);

            % parse the inputs
            parse(p,varargin{:});
            % store the inputs as a struct
            inp = p.Results;
            % get the fieldnames of the edge functions
            inp = orderfields(inp,{'expand','velocity_model'});
            flds = fieldnames(inp);
            for i = 1 : numel(flds)
                type = flds{i};
                switch type
                    case('expand')
                        obj.expand = inp.(flds{i}); 
                        if ~isempty(obj.expand) 
                            obj.expand = inp.(flds{i}); 
                        end
                    case('velocity_model')
                        obj.fname = inp.(flds{i});
                        if ~isempty(obj.fname)
                            obj.fname = inp.(flds{i});
                            obj = ReadVelocityData(obj);
                        end
                end
            end
        end
        
        % getters
        function F=GetFvp(obj), F = obj.Fvp; end
        
        function ny=GetNy(obj), ny = obj.ny; end
        
        function ny=GetNx(obj), ny = obj.ny; end
        
        function nz=GetNz(obj), nz = obj.nz; end
        
        function dim=GetDim(obj), dim = obj.dim; end
        
        function gsp=GetGridspace(obj), gsp = obj.gridspace; end

        % plotting
        function plot(obj)
            if(obj.GetDim == 2)
                [yg,zg]=CreateStructGrid(obj) ;
                tmp=obj.Fvp(yg,zg);
                skip=1 ; % save memory and time by skipping
                figure;
                pcolor(zg(1:skip:end,1:skip:end)*obj.gridspace,...
                    yg(1:skip:end,1:skip:end)*obj.gridspace,...
                    tmp(1:skip:end,1:skip:end)) ;
                shading interp;
                set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
                xlabel('X-position (m)');
                ylabel('Z-position/depth (m)');
                cb=colorbar; ylabel(cb,'P-wave speed (m/s)') ;
                set(gca,'FontSize',16) ;
                caxis([1000 5000])
            elseif obj.GetDim == 3
                [xg,yg,zg]=CreateStructGrid3D(obj);
                tmp=obj.Fvp(xg,yg,zg);
                figure; scatter3(xg(:),yg(:),zg(:),5,tmp(:) ) ;
                set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
                xlabel('X-position (m)');
                ylabel('Y-position (m)');
                zlabel('Z-position/depth (m)');
                cb=colorbar; ylabel(cb,'P-wave speed (m/s)') ;
                set(gca,'FontSize',16) ;
                axis equal; axis tight;
            end
        end
        
        function [yg,zg]=CreateStructGrid(obj)
            [yg,zg] = ndgrid(obj.x0y0(1) + (0:obj.ny-1)'*obj.gridspace, ...
                obj.x0y0(2) + (0:obj.nz-1)'*obj.gridspace);
        end
                
        function [xg,yg,zg]=CreateStructGrid3D(obj)
            [xg,yg,zg] = ndgrid((0:obj.nx-1)'*obj.gridspace,...
                obj.x0y0(1) + (0:obj.ny-1)'*obj.gridspace, ...
                obj.x0y0(2) + (0:obj.nz-1)'*obj.gridspace);
        end
        
        function [vpOnNodes]=InterpVP(obj,P)
            f=GetFvp(obj); 
            vpOnNodes=f(P(:,1:2)); 
        end
        
        
    end % end non-static public methods
    
    methods(Access=private)
        
        function obj = ReadVelocityData(obj)
            if exist(obj.fname, 'file') == 2
                obj.dim=ncread(obj.fname,'dim');
                obj.gridspace=ncread(obj.fname,'gridspace') ;
                % File exists.
                if obj.dim == 2
                    tmp = ncread(obj.fname,'vp');
                    
                    % determine size of grid
                    obj.ny = ncread(obj.fname,'ny');
                    obj.nz = ncread(obj.fname,'nz');
                    
                    % determine bottom front corner coordinate
                    obj.x0y0(1:obj.dim) = ncread(obj.fname,'x0y0z0');
                    
                    [zg,yg]=CreateStructGrid(obj);
                    
                    obj.bbox = [min(zg(:)) max(zg(:))
                        min(yg(:)) max(yg(:))];
                    
                    obj.Fvp=griddedInterpolant(zg,yg,tmp,'linear','none') ;

                    % if expansion is enabled, we must also expand the
                    % velocity model with a constant velocity 
                    if obj.expand~=0
                        
                        add_ny_ = ceil(obj.expand/obj.gridspace);
                        add_nz_ = add_ny_ ;
                        
                        y_exp = add_ny_ * obj.gridspace;
                        
                        x0y0_(1) = obj.x0y0(1)  ;
                        x0y0_(2) = obj.x0y0(2) - y_exp;
                        
                        nz_ = obj.nz + 2*add_nz_ ;
                        ny_ = obj.ny + add_ny_ ;
                        
                        [zg_,yg_] = ndgrid(x0y0_(1) + (0:ny_-1)'*obj.gridspace, ...
                            x0y0_(2) + (0:nz_-1)'*obj.gridspace);
                        
                        bbox_ = [min(zg_(:)) max(zg_(:))
                            min(yg_(:)) max(yg_(:))];
                        
                        % interpolant the smaller velocity model onto the
                        % expanded domain 
                        vp_ = obj.Fvp(zg_,yg_);
                        vp_(isnan(vp_))= NaN ; % FLAG VALUE 
                        obj.Fvp=griddedInterpolant(zg_,yg_,vp_) ;

                        obj.bbox = bbox_ ; 
                        obj.ny = ny_ ; 
                        obj.nz = nz_ ; 
                        obj.x0y0 = x0y0_; 
                        
                        % also have the code output a new Vp_EXACT Vp_GUESS
                        % to run the FWI code 
                    end
                    
                    clearvars yg zg tmp;
                    
                    disp(['INFO: SUCCESFULLY READ IN FILE ',obj.fname]);
                else
                    
                    tmp = ncread(obj.fname,'vp');
                    
                    % determine size of grid
                    obj.nx = ncread(obj.fname,'nx');
                    obj.ny = ncread(obj.fname,'ny');
                    obj.nz = ncread(obj.fname,'nz');

                    % determine bottom front corner coordinate
                    obj.x0y0(1:obj.dim) = ncread(obj.fname,'x0y0z0');
                    
                    [xg,yg,zg]=CreateStructGrid3D(obj);
                    
                    obj.bbox = [min(xg(:)) max(xg(:))
                        min(yg(:)) max(yg(:))
                        min(zg(:)) max(zg(:))];
                    
                    obj.Fvp=griddedInterpolant(xg,yg,zg,tmp) ;
                    
                    clearvars xg yg zg tmp;
                    
                    disp(['INFO: SUCCESFULLY READ IN FILE ',obj.fname]);
                end
            else
                % File does not exist.
                error(['ERROR: FILE CANNOT BE LOCATED ',obj.fname]);
            end
        end
        
    end
    
end

