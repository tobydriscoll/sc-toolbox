classdef  (InferiorClasses = {?double}) diskmap < scmap
    
    properties
        prevertex = [];
        constant = [];
        qdata = [];
        accuracy = [];
        center = [];
    end
    
    methods
        function map = diskmap(varargin)
            %DISKMAP Schwarz-Christoffel disk map object.
            %   DISKMAP(P) constructs a Schwarz-Christoffel disk map object for the
            %   polygon P. The parameter problem is solved using default options for
            %   the prevertices and the multiplicative constant.
            %
            %   DISKMAP(P,OPTIONS) uses an options structure of the type created by
            %   SCMAPOPT in solving the parameter problem.
            %
            %   DISKMAP(P,Z) creates a diskmap object having the given prevertices Z
            %   (the mulitiplicative constant is found automatically). There is no
            %   checking to ensure that the prevertices are consistent with the
            %   given polygon. DISKMAP(P,Z,C) also uses the supplied constant. An
            %   OPTIONS argument can be added, although only the error tolerance
            %   will be used.
            %
            %   DISKMAP(M), where M is a diskmap object, just returns M.
            %
            %   DISKMAP(M,P) returns a new diskmap object for the polygon P using
            %   the options in diskmap M. The prevertices of M will be used as the
            %   starting guess for the parameter problem of the new map. Thus P
            %   should properly be a perturbation of the polygon for M. An OPTIONS
            %   structure may be added to override options in M.
            %
            %   DISKMAP(Z,ALPHA) creates a map using the given prevertices and the
            %   interior polygon angles described by ALPHA (see POLYGON help). The
            %   image polygon is deduced by computing S-C integrals assuming a
            %   multiplicative constant of 1. DISKMAP(Z,ALPHA,C) uses the given
            %   constant instead.
            %
            %   See also SCMAPOPT, classes POLYGON, SCMAP.
            
            %   Copyright 1998--2001 by Toby Driscoll.
            
            
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
                case 'diskmap'
                    oldmap = varargin{1};
                    % Continuation of given map to given polygon
                    poly = varargin{2};
                    opt = scmapopt(options(oldmap));
                    z0 = prevertex(oldmap);
                    if length(z0) ~= length(poly)
                        msg = 'Polygon %s must have the same length as that in %s.';
                        error(sprintf(msg,inputname(2),inputname(1)))
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
                            z = z(:);
                        elseif length(arg) == 1
                            c = arg;
                        else
                            msg = 'Unable to parse argument ''%s''.';
                            error(sprintf(msg,inputname(j)))
                        end
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
                            error(sprintf(msg,inputname(j)))
                        end
                    end
                    
                otherwise
                    msg = 'Expected ''%s'' to be a polygon, a diskmap, or prevertices.';
                    error(sprintf(msg,inputname(1)))
            end % input switch
            
            % Retrieve options
            opt = scmapopt(opt);
            
            % Take actions based on what needs to be filled in
            
            if isempty(z)
                % Solve parameter problem
                % Apply SCFIX to enforce solver rules
                [w,beta] = scfix('d',vertex(poly),angle(poly)-1);
                poly = polygon(w,beta+1);             % in case polygon was altered
                
                [z,c,qdata] = dparam(w,beta,opt.InitialGuess,opt);
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
                I = dquad(z(1),mid,1,z,beta,qdata) - dquad(z(idx),mid,idx,z,beta,qdata);
                c = diff(w([1 idx]))/I;
            end
            
            map = map@scmap(poly,opt);
            
            map.prevertex = z;
            map.constant = c;
            map.qdata = qdata;
            
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
    
    methods (Static)
        % needed by CR inverse mapping
        function out = imapfun(varargin)
            out = dimapfun(varargin{:});
        end
        
        % Needed in @hplmap/private/hp2disk
        function out = dquad(varargin)
            out = dquad(varargin{:});
        end
    end
        
end
