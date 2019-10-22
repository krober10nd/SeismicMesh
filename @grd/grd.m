classdef grd
    % A container for a mesh
    
    properties(Access=private)
        p % points
        t % facets
    end
    
    methods(Access=public)
        function obj = grd(p,t)
            if  nargin > 1
              % Constructor given p and t 
              obj.p = p; 
              obj.t = t; 
            else
               obj.p = []; 
               obj.t = [];
            end
        end
        
        function p=GetP(obj), p=obj.p; end
        
        function t=GetT(obj); t=obj.t; end
    end
    
    
end

