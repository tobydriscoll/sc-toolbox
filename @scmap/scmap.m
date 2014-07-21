classdef (InferiorClasses = {?double}) scmap
    
    properties
        polygon = []
        options = []
    end
    
    methods
        function map = scmap(poly,opt)
            %SCMAP Construct generic Schwarz-Christoffel map object.
            %   SCMAP(P) creates a generic parent scmap object whose target polygon
            %   is given by P. SCMAP(P,OPTIONS) accepts an options structure
            %   produced by SCMAPOPT.
            %
            %   SCMAP(M), where M is already an scmap object, returns M. SCMAP by
            %   itself returns an object with empty polygon and default options.
            %
            %   You do not need to create an scmap object directly. Use one of the
            %   specific child classes instead.
            %
            %   See also SCMAPOPT, and classes DISKMAP, HPLMAP, EXTERMAP, STRIPMAP,
            %   RECTMAP, CRDISKMAP, CRRECTMAP.
            
            %   Copyright 1998 by Toby Driscoll.
                        
            if nargin == 0
                % Leave everything empty
                poly = [];
                opt = [];
            else
                
                % Branch based on class of first argument
                switch class(poly)
                    case 'polygon'
                        if nargin == 1
                            opt = [];
                        end
                    otherwise
                        msg = 'Expected ''%s'' to be of class polygon or scmap.';
                        error(sprintf(msg,inputname(1)))
                end
                
            end
            
            map.polygon = poly;
            map.options = sctool.scmapopt(opt);
                        
        end
    end
    
end