classdef scmapdiff
    
    properties
        themap
    end
    
    methods
        
        function df = scmapdiff(f)
            %SCMAPDIFF Derivative of a Schwarz-Christoffel map.
            %   SCMAPDIFF(F) returns a dummy object that represents the derivative
            %   of the SC map F. The only thing possible with this object is to EVAL
            %   it.
            %
            %   See also SCMAPDIFF/EVAL, SCMAPDIFF/SUBSREF.
            
            %   Copyright 1998-2001 by Toby Driscoll.
            %   $Id: scmapdiff.m 160 2001-07-20 14:28:22Z driscoll $
            
            if nargin > 0
                df.themap = f;
            else
                df.themap = [];
            end
            
        end
    end
end