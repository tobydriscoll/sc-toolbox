classdef  (InferiorClasses = {?double}) extermap < scmap
    
    properties
        prevertex = [];
        constant = [];
        qdata = [];
        accuracy = [];
    end
    
    methods
        
        function map = extermap(varargin)
            %EXTERMAP Schwarz-Christoffel exterior map object.
            %   EXTERMAP(P) constructs a Schwarz-Christoffel exterior map object for
            %   the polygon P. The parameter problem is solved using default options
            %   for the prevertices and the multiplicative constant.
            %
            %   EXTERMAP(P,OPTIONS) uses an options structure of the type created by
            %   SCMAPOPT in solving the parameter problem.
            %
            %   EXTERMAP(P,Z) creates a extermap object having the given prevertices
            %   Z (the mulitiplicative constant is found automatically).
            %   EXTERMAP(P,Z,C) also uses the given constant. An OPTIONS argument
            %   can be added, although only the error tolerance will be used.
            %
            %   EXTERMAP(M), where M is a extermap object, just returns M.
            %
            %   EXTERMAP(M,P) returns a new extermap object for the polygon P using
            %   the options in extermap M. The prevertices of M will be used as the
            %   starting guess for the parameter problem of the new map. Thus P
            %   should properly be a perturbation of the polygon for M. An OPTIONS
            %   structure may also be given to override options in M.
            %
            %   EXTERMAP(Z,ALPHA) creates a map using the given prevertices and the
            %   interior polygon angles described by ALPHA (see POLYGON help). The
            %   image polygon is deduced by computing S-C integrals assuming a
            %   multiplicative constant of 1. EXTERMAP(Z,ALPHA,C) uses the given
            %   constant instead. Note that not every pairing of prevertices and
            %   angles produces a single-valued map; you must have SUM((ALPHA-1)./Z)
            %   equal to zero. Also, Z is given counterclockwise around the unit
            %   circle, but ALPHA should be clockwise with respect to the interior
            %   of the polygon.
            %
            %   See also SCMAPOPT, classes POLYGON, SCMAP.
            
            %   Copyright 1998-2001 by Toby Driscoll.
            %   $Id: extermap.m 129 2001-05-07 15:04:13Z driscoll $
            
            % Initialize with empties
            poly = [];
            alpha = [];
            z = [];
            c = [];
            opt = [];
            qdata = [];
            
            import sctool.*
            
            % Branch based on class of first argument
            switch class(varargin{1})
                case 'extermap'
                    oldmap = varargin{1};
                    % Continuation of given map to given polygon
                    poly = varargin{2};
                    opt = scmapopt(options(oldmap));
                    z0 = oldmap.prevertex;
                    if length(z0) ~= length(poly)
                        msg = 'Polygon %s must have the same length as that in %s.';
                        error(msg,inputname(2),inputname(1))
                    end
                    if nargin > 2
                        opt = scmapopt(opt,varargin{3});
                    end
                    opt = scmapopt(opt,'initial',z0);
                    
                case 'polygon'
                    poly = varargin{1};
                    % Parse optional arguments
                    for j = 2:length(varargin)
                        arg = varargin{j};
                        % Each arg is an options struct, z, or c
                        if isa(arg,'struct')
                            opt = arg;
                        elseif length(arg) == length(poly)
                            z = arg;
                            % We will have to flip vertices to get correct orientation
                            z = flipud(z(:));
                        elseif isscalar(arg)
                            c = arg;
                        else
                            msg = 'Unable to parse argument ''%s''.';
                            error(msg,inputname(j+1))
                        end
                    end
                    
                case 'double'
                    % Args are the prevertex vector, then angle vector
                    z = varargin{1}(:);
                    alpha = varargin{2};
                    poly = polygon(NaN*alpha*1i,alpha);
                    c = 1;
                    % Check residue of integrand to see if compatible
                    if abs(sum((alpha-1)./z)) > 1e-8
                        error('Map is not single-valued')
                    end
                    for j = 3:length(varargin)
                        if isa(varargin{j},'struct')
                            opt = varargin{j};
                        elseif isscalar(varargin{j})
                            c = varargin{j};
                        else
                            msg = 'Unable to parse argument ''%s''.';
                            error(msg,inputname(j+1))
                        end
                    end
                    
                otherwise
                    msg = 'Expected ''%s'' to be of class polygon or extermap.';
                    error(msg,inputname(1))
            end % switch
            
            
            % Retrieve options
            opt = scmapopt(opt);
            
            % Take actions based on what needs to be filled in
            
            if isempty(z)
                % Find prevertices
                % Apply SCFIX to enforce solver rules
                w = flipud(vertex(poly));
                beta = 1 - flipud(angle(poly));
                [w,beta] = scfix('de',w,beta);
                poly = polygon(flipud(w),1-flipud(beta));
                
                z0 = opt.InitialGuess;
                tol = opt.Tolerance;
                [z,c,qdata] = deparam(w,beta,z0,opt);
            end
            
            if isempty(qdata)
                % Base quadrature accuracy on given options
                nqpts = ceil(-log10(opt.Tolerance));
                beta = 1 - flipud(angle(poly));
                qdata = scqdata(beta,nqpts);
            end
            
            if isempty(c)
                % Find constant
                w = flipud(vertex(poly));
                beta = 1 - flipud(angle(poly));
                mid = z(1)*exp(1i*angle(z(2)/z(1))/2);
                I = dequad(z(1),mid,1,z,beta,qdata) - dequad(z(2),mid,2,z,beta,qdata);
                c = diff(w(1:2))/I;
            end
            
            map = map@scmap(poly,opt);
            
            map.prevertex = z;
            map.constant = c;
            map.qdata = qdata;
 
            % If the polygon was not known, find it from the map
            if any(isnan(vertex(poly)))
                poly = forwardpoly(map);
                map = extermap(poly,z,opt);
            end
            
            % Now fill in apparent accuracy
            map.accuracy = accuracy(map);
            
        end
    end
end
