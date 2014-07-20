classdef  (InferiorClasses = {?double}) stripmap < scmap
    
    properties
        prevertex = []
        constant = []
        qdata = []
        accuracy = []
    end
    
    methods
        
        function map = stripmap(poly,varargin)
            %STRIPMAP Schwarz-Christoffel strip map object.
            %   STRIPMAP(P,ENDIDX) constructs a Schwarz-Christoffel strip map object
            %   for the polygon P. ENDIDX is a two-vector containing the indices of
            %   the vertices that are the images of the left and right ends of the
            %   strip. STRIPMAP(P) requires you to choose these vertices
            %   graphically. The parameter problem is solved using default options
            %   for the prevertices and the multiplicative constant.
            %
            %   STRIPMAP(P,ENDIDX,OPTIONS) or STRIPMAP(P,OPTIONS) uses an options
            %   structure of the type created by SCMAPOPT in solving the parameter
            %   problem.
            %
            %   STRIPMAP(P,Z) creates a stripmap object having the given prevertices
            %   Z (the mulitiplicative constant is found automatically). Z must
            %   include -Inf and Inf, the ends of the strip.  STRIPMAP(P,Z,C) also
            %   uses the given constant. An OPTIONS argument can be added, although
            %   only the error tolerance will be used.
            %
            %   STRIPMAP(M), where M is a stripmap object, just returns M.
            %
            %   STRIPMAP(M,P) returns a new stripmap object for the polygon P using
            %   the options in stripmap M. The prevertices of M will be used as the
            %   starting guess for the parameter problem of the new map. Thus P
            %   should properly be a perturbation (continuation) of the polygon for
            %   M. An OPTIONS structure may be added to override options in M. There
            %   is no opportunity to change the end indices.
            %
            %   See also SCMAPOPT, and classes POLYGON, SCMAP.
            
            %   Copyright 1998 by Toby Driscoll.
            
            % Assign empties to optional args
            endidx = [];
            z = [];
            c = [];
            opt = [];
            
            % Branch based on class of first argument
            switch class(poly)
                
                case 'stripmap'
                    oldmap = poly;
                    % Continuation of given map to given polygon
                    poly = varargin{1};
                    z0 = oldmap.prevertex;
                    if length(z0) ~= length(poly)
                        msg = 'Polygon %s must have the same length as that in %s.';
                        error(msg,inputname(2),inputname(1))
                    end
                    opt = scmapopt(oldmap.scmap);
                    if nargin > 2
                        opt = scmapopt(opt,varargin{2});
                    end
                    opt = scmapopt(opt,'initial',z0);
                    endidx(1) = find( isinf(z0) & (z0 < 0) );
                    endidx(2) = find( isinf(z0) & (z0 > 0) );
                    
                case 'polygon'
                    % Parse optional arguments
                    for j = 1:length(varargin)
                        arg = varargin{j};
                        % Each arg is the end specifier, an options struct, z, or c
                        if isa(arg,'struct')
                            opt = arg;
                        elseif (length(arg) == 2) & all(round(arg) == arg)
                            endidx = arg;
                        elseif length(arg) == length(poly)
                            z = arg;
                            z = z(:);
                        elseif length(arg) == 1
                            c = arg;
                        else
                            msg = 'Unable to parse argument ''%s''.';
                            error(msg,inputname(j+1))
                        end
                    end
                    
                otherwise
                    msg = 'Expected ''%s'' to be of class polygon or stripmap.';
                    error(msg,inputname(1))
                    
            end % switch
            
            % Retrieve options
            opt = scmapopt(opt);
            
            
            % Get data for the low-level functions
            w = vertex(poly);
            n = length(w);
            beta = angle(poly) - 1;
            
            % Request endidx
            if isempty(endidx)
                msg = 'Select the vertices that map to the ends of the strip.';
                endidx = scselect(w,beta,2,'Select ends',msg);
            end
            
            % Find prevertices if necessary
            if isempty(z)
                
                % Apply SCFIX to enforce solver rules
                [w,beta,endidx] = scfix('st',w,beta,endidx);
                poly = polygon(w,beta+1);
                
                [z,c,qdata] = stparam(w,beta,endidx,opt.InitialGuess,opt);
                
            else
                
                atinf = isinf(z);
                if (sum(atinf & (z < 0)) ~= 1) | (sum(atinf & (z > 0)) ~= 1)
                    error('Supplied prevertices must include -Inf and Inf')
                end
                % Base quadrature accuracy on given options
                nqpts = ceil(-log10(opt.Tolerance));
                qdata = scqdata(beta,nqpts);
                
            end
            
            % Find constant if necessary
            if isempty(c)
                k = find(isinf(z) & (z < 0));
                mid = mean(z(k+1:k+2));
                I = stquad(z(k+1),mid,k+1,z,beta,qdata) - ...
                    stquad(z(k+2),mid,k+2,z,beta,qdata);
                c = diff(w(k+1:k+2))/I;
            end
            
            map = map@scmap(poly,opt);
            map.prevertex = z;
            map.constant = c;
            map.qdata = qdata;
            
            % Now fill in apparent accuracy
            map.accuracy = accuracy(map);
            
        end
    end
end

