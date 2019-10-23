classdef geodata
    %GEODATA Read in a velocity model from a SEGY file
    
    %-----------------------------------------------------------
    %   Keith Roberts   : 2019 --
    %   Email           : krober@usp.br
    %   Last updated    : 10/20/2019
    %-----------------------------------------------------------
    %
    properties(Access=private)
        fname % path to segy seismic data.
        x0y0  % top left corner coordinate
        Fvp   % gridded interpolant of seismic wavespeed in km/s
        ny    % number of grid points in y direction
        nz    % number of grid points in z direction
        gridspace    % grid space (in m)
    end
    
    properties(Access=public)
        bbox % domain corners left, right, bottom, top.
    end
        
    methods(Access=public)
        % default class constructor
        % GEODATA construct the default class.
        % loads in a 2D velocity model from a segy file.
        function obj = geodata(varargin)
            %
            p = inputParser;
            
            defval=0;
            
            % add name/value pairs
            addOptional(p,'segy',defval);
            addOptional(p,'gridspace',defval);

            % parse the inputs
            parse(p,varargin{:});
            % store the inputs as a struct
            inp = p.Results;
            % get the fieldnames of the edge functions
            inp = orderfields(inp,{'gridspace','segy'});
            flds = fieldnames(inp);
            for i = 1 : numel(flds)
                type = flds{i};
                switch type
                    case('segy')
                        obj.fname = inp.(flds{i});
                        if ~isempty(obj.fname)
                            obj.fname = inp.(flds{i});
                            obj = ReadVelocityData(obj);
                        end
                    case('gridspace')
                        obj.gridspace = inp.(flds{i}); 
                        if obj.gridspace~=0
                            obj.gridspace = inp.(flds{i});
                        else
                            error('Please pass a gridspace in meters to geodata!');
                        end
                end
            end
        end
        
        % getters
        function F=GetFvp(obj), F = obj.Fvp; end
        
        function ny=GetNy(obj), ny = obj.ny; end
        
        function nz=GetNz(obj), nz = obj.nz; end

        function gsp=GetGridspace(obj), gsp = obj.gridspace; end

        % plotting
        function [axH]=plot(obj)
            [yg,zg]=CreateStructGrid(obj) ;
            tmp=obj.Fvp(yg,zg);
            skip=5 ; % save memory and time by skipping
            figure;
            axH=pcolor(yg(1:skip:end,1:skip:end)*obj.gridspace,...
                zg(1:skip:end,1:skip:end)*obj.gridspace,...
                tmp(1:skip:end,1:skip:end)) ;
            shading interp;
            set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
            xlabel('Y-position (m)');
            ylabel('Z-position/depth (m)');
            cb=colorbar; ylabel(cb,'P-wave speed (km/s)') ;
            set(gca,'FontSize',16) ;
        end
        
        function [xg,yg]=CreateStructGrid(obj)
            [xg,yg] = ndgrid(obj.x0y0(1) + (0:obj.ny-1)'*obj.gridspace, ...
                obj.x0y0(2) + (0:obj.nz-1)'*obj.gridspace);
        end
        
    end % end non-static public methods
    
    methods(Access=private)
        
        function obj = ReadVelocityData(obj)
            if exist(obj.fname, 'file') == 2
                % File exists.
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
                % File does not exist.
                error(['ERROR: FILE CANNOT BE LOCATED ',segy]);
            end
        end
        
    end
    
end

