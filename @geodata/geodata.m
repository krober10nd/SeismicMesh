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

            % parse the inputs
            parse(p,varargin{:});
            % store the inputs as a struct
            inp = p.Results;
            % get the fieldnames of the edge functions
            inp = orderfields(inp,{'velocity_model'});
            flds = fieldnames(inp);
            for i = 1 : numel(flds)
                type = flds{i};
                switch type
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
            [yg,zg]=CreateStructGrid(obj) ;
            tmp=obj.Fvp(yg,zg);
            skip=5 ; % save memory and time by skipping
            figure;
            pcolor(yg(1:skip:end,1:skip:end)*obj.gridspace,...
                zg(1:skip:end,1:skip:end)*obj.gridspace,...
                tmp(1:skip:end,1:skip:end)) ;
            shading interp;
            set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
            xlabel('Y-position (m)');
            ylabel('Z-position/depth (m)');
            cb=colorbar; ylabel(cb,'P-wave speed (m/s)') ;
            set(gca,'FontSize',16) ;
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
        
        
    end % end non-static public methods
    
    methods(Access=private)
        
        function obj = ReadVelocityData(obj)
            if exist(obj.fname, 'file') == 2
                obj.dim=ncread(obj.fname,'dim'); 
                obj.gridspace=ncread(obj.fname,'gridspace') ; 
                % File exists.
                if obj.dim == 2
                    tmp=ReadSegy(obj.fname)';
                    [obj.ny,obj.nz]=size(tmp) ;
                    obj.x0y0=[0,0];
                    [yg,zg]=CreateStructGrid(obj);
                    obj.bbox = [min(yg(:)) max(yg(:))
                        min(zg(:)) max(zg(:))];
                    obj.Fvp=griddedInterpolant(yg,zg,tmp) ;
                    clearvars yg zg tmp;
                    disp(['INFO: SUCCESFULLY READ IN FILE',obj.fname]);
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

