classdef (InferiorClasses = {?double}) riesurfmap < scmap
    
    properties
        branch = [];
        prevertex = [];
        prebranch = [];
        constant = [];
        qdata = [];
        accuracy = [];
        center = [];
    end
    
    methods
        function map = riesurfmap(varargin)
            %RIESURFMAP Schwarz-Christoffel map to Riemann surface.
            
            %   RIESURFMAP(P,BRANCH) constructs a Schwarz-Christoffel map to the
            %   region bounded by polygon P having one or more branch points in
            %   BRANCH. One must have
            %
            %      sum( angle(P)-1 ) == -2*(1+length(BRANCH)),
            %
            %   and the branch points must be on multiply covered parts of the
            %   plane. No effort is made to check these conditions. The parameter
            %   problem is solved for the prevertices and the multiplicative
            %   constant using default options.
            %
            %   RIESURFMAP(P,BRANCH,OPTIONS) uses an options structure of the type
            %   created by SCMAPOPT in solving the parameter problem.
            %
            %   RIESURFMAP(P,BRANCH,Z,ZB) creates a map having the given prevertices
            %   Z (the mulitiplicative constant is found automatically) and
            %   prebranch points ZB. There is no checking to ensure that these data
            %   are consistent with the given region. RIESURFMAP(P,BRANCH,Z,ZB,C)
            %   also uses the supplied constant. An OPTIONS argument can be added,
            %   although only the error tolerance will be used.
            %
            %   RIESURFMAP(F), where F is a riesurfmap, just returns F.
            %
            %   RIESURFMAP(F,P,BRANCH) returns a new map for the polygon P using the
            %   options in riesurfmap F. The prevertices of F will be used as the
            %   starting guess for the parameter problem of the new map. Thus P
            %   should properly be a perturbation of the polygon for F. An OPTIONS
            %   structure may be added to override options in F.
            %
            %   RIESURFMAP(Z,ZB,ALPHA) creates a map using the given
            %   prevertices/prebranches and the interior polygon angles described by
            %   ALPHA (see POLYGON help). The image polygon is deduced by computing
            %   S-C integrals, assuming a multiplicative constant of unity.
            %   RIESURFMAP(Z,ALPHA,C) uses the given constant instead.
            %
            %   See also SCMAPOPT, classes POLYGON, SCMAP.
            
            %   Copyright 2002 by Toby Driscoll.
            
            % Initialize with empties
            poly = [];
            branch = [];
            alpha = [];
            z = [];
            zb = [];
            c = [];
            opt = [];
            qdata = [];
            import sctool.*
            
            % Branch based on class of first argument
            switch class(varargin{1})
                case 'riesurfmap'
                    oldmap = varargin{1};
                    % Continuation of given map to given polygon
                    poly = varargin{2};
                    branch = varargin{3};
                    opt = scmapopt(oldmap);
                    z0 = prevertex(oldmap);
                    zb = prebranch(oldmap);
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
                    branch = varargin{2};
                    % Parse optional arguments
                    j = 3;
                    while j <= length(varargin)
                        arg = varargin{j};
                        % Is it an options struct, z, or c?
                        if isa(arg,'struct')
                            opt = arg;
                        elseif length(arg) == length(poly)
                            z = arg;
                            z = z(:);
                            zb = varargin{j+1};
                            j = j+1;
                        elseif length(arg) == 1
                            c = arg;
                        else
                            msg = 'Unable to parse argument ''%s''.';
                            error(msg,inputname(j+1))
                        end
                        j = j+1;
                    end
                    
                case 'double'
                    % Args are the prevertex vector, then angle vector
                    z = varargin{1};
                    alpha = varargin{2};
                    poly = polygon(NaN*alpha*1i,alpha);
                    c = 1;
                    for j = 3:length(varargin)
                        if isa(varargin{j},'struct')
                            opt = varargin{j};
                        elseif length(varargin{j})==1
                            c = varargin{j};
                        else
                            msg = 'Unable to parse argument ''%s''.';
                            error(msg,inputname(j+1))
                        end
                    end
                    
                otherwise
                    msg = 'Expected ''%s'' to be a polygon, a riesurfmap, or prevertices.';
                    error(msg,inputname(1))
            end % input switch
            
            % Retrieve options
            opt = scmapopt(opt);
            
            % Take actions based on what needs to be filled in
            
            if isempty(z)
                % Solve parameter problem
                % Apply SCFIX to enforce solver rules
                [w,beta] = scfix('d',vertex(poly),angle(poly)-1);
                poly = polygon(w,beta+1);             % in case polygon was altered
                
                [z,zb,c,qdata] = rsparam(w,beta,branch,opt.InitialGuess,opt);
            end
            
            if isempty(qdata)
                % Base accuracy of quadrature on given options
                nqpts = ceil(-log10(opt.Tolerance));
                qdata = scqdata(angle(poly)-1,nqpts);
            end
            
            if isempty(c)
                % Find constant
                w = vertex(poly);
                beta = angle(poly)-1;
                idx = 1 + find(~isinf(z(2:end)), 1 );
                mid = mean(z([1 idx]));
                I = rsquad(z(1),z(2),1,2,z,beta,zb,qdata);
                c = diff(w([1 idx]))/I;
            end
            
            map = map@scmap(poly,opt);
            
            map.branch = branch;
            map.prevertex = z;
            map.prebranch = zb;
            map.constant = c;
            map.qdata = qdata;
            xs           
            % If the polygon was not known, find it from the map
            if any(isnan(vertex(poly)))
                poly = forwardpoly(map);
                map.scmap = scmap(poly,opt);
            end
            
            % Find conformal center
            map.center = center(map);
            
            % Fill in apparent accuracy
            map.accuracy = accuracy(map);
            
        end
        
    end
    
end