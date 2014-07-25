classdef annulusmap
    %ANNULUSMAP Conformal map to doubly connected polygonal region.
    %   ANNULUSMAP(POLYOUTER,POLYINNER) creates a Schwarz-Christoffel annulus
    %   map to the region bounded by polygons POLYOUTER and POLYINNER.
    %
    %   ANNULUSMAP(POLYOUTER,POLYINNER,'truncate') will truncate any unbounded
    %   sides. Otherwise an error is thrown for unbounded regions.
    
    %   Copyright by Toby Driscoll, 2014
    %   Written by Alfa Heryudono, 2003 and Toby Driscoll, 2014.
    
    properties
        boundary
        M, N, Z0, Z1, ALFA0, ALFA1
        isUnbounded
        qwork
        u, c, w0, w1, phi0, phi1
    end
    
    properties (Dependent)
        % backward compatibility
        ISHAPE
    end
    
    methods
        
        function map = annulusmap(outerPolygon,innerPolygon,varargin)
            
            map.boundary = { outerPolygon, innerPolygon };
            
            if isinf(outerPolygon)
                if (nargin < 3) || ~isequal(varargin{1},'truncate')
                    error('Region must be bounded, or extra option given.')
                else
                    outerPolygon = truncate(outerPolygon);
                    map.isUnbounded = true;
                end
            else
                map.isUnbounded = false;
            end
            
            map.M = length(outerPolygon); map.N = length(innerPolygon);
            map.Z0 = vertex(outerPolygon).';
            map.Z1 = vertex(innerPolygon).';
            map.ALFA0 = angle(outerPolygon)'; 
            map.ALFA1 = 2 - angle(innerPolygon)';
            
            % Generate the Gauss-Jacobi weights & nodes
            nptq = 8;
            map.qwork = qinit(map,nptq);
            
            % Solve parameter problem.
            iguess = 0;  %(0 nonequally spaced guess or 1 equally spaced guess)
            linearc = 1; % (0 line path or 1 circular path)
            [map.u,map.c,map.w0,map.w1,map.phi0,map.phi1] = ...
                dscsolv(iguess,nptq,map.qwork,map.isUnbounded,linearc,map);
            
        end
        
        function I = get.ISHAPE(map)
            I = map.isUnbounded;
        end
    end
end
