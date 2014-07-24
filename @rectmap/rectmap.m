classdef  (InferiorClasses = {?double}) rectmap < scmap
    
    properties
        prevertex = []
        constant = []
        qdata = []
        stripL = []
        accuracy = []
    end
    
    methods
        
        function map = rectmap(poly,varargin)
            %RECTMAP Schwarz-Christoffel rectangle map object.
            %   RECTMAP(P,CORNER) constructs a Schwarz-Christoffel rectangle map
            %   object for the polygon P. CORNER is a four-vector containing the
            %   indices of the vertices that are the images of the rectangle's
            %   corners. These indices should be specified in counterclockwise
            %   order, and the first two should be the endpoints of a long side of
            %   the rectangle.
            %
            %   RECTMAP(P) requires you to choose the corners graphically.
            %
            %   RECTMAP(P,CORNER,OPTIONS) or RECTMAP(P,OPTIONS) uses an options
            %   structure of the type created by SCMAPOPT in solving the parameter
            %   problem.
            %
            %   RECTMAP(P,Z,C,L) constructs a rectmap using the explicitly given
            %   prevertices Z, constant C, and strip length L.
            %
            %   RECTMAP(M), where M is a rectmap object, just returns M.
            %
            %   RECTMAP(M,P) returns a new rectmap object for the polygon P using
            %   the options in rectmap M. The prevertices of M will be used as the
            %   starting guess for the parameter problem of the new map. Thus P
            %   should properly be a perturbation (continuation) of the polygon for
            %   M. An OPTIONS structure may be added to override options in M. There
            %   is no opportunity to change the corner indices.
            %
            %   See also SCMAPOPT.
            
            %   Copyright 1998--2001 by Toby Driscoll.
         
            % Assign empties to optional args
            corner = [];
            z = [];
            c = [];
            L = [];
            opt = [];
            import sctool.*
            
            % Branch based on class of first argument
            switch class(poly)
                
                case 'rectmap'
                    oldmap = poly;
                    % Continuation of given map to given polygon
                    poly = varargin{1};
                    % Use strip prevertices as initial guess
                    %z0 = M.strip_prevertex;
                    zr = oldmap.prevertex;
                    z0 = r2strip(zr,zr,oldmap.stripL);
                    if length(z0) ~= length(poly)
                        msg = 'Polygon %s must have the same length as that in %s.';
                        error(msg,inputname(2),inputname(1))
                    end
                    opt = scmapopt(oldmap);
                    if nargin > 2
                        opt = scmapopt(opt,varargin{2});
                    end
                    opt = scmapopt(opt,'initial',z0);
                    corner = corners(oldmap);
                    
                case 'polygon'
                    % Parse optional arguments
                    for j = 1:length(varargin)
                        arg = varargin{j};
                        % Each arg is the corner vector, an options struct, or z/c/L
                        if isa(arg,'struct')
                            opt = arg;
                        elseif (length(arg) == 4) && all(round(arg) == arg)
                            corner = arg;
                        elseif length(arg) == length(poly)
                            % In this case, immediately pick up the two constants as well
                            z = arg;
                            z = z(:);
                            if j~=1 or length(varargin)~=3
                                error('Incorrectly specified prevertices, c, and L.')
                            else
                                c = varargin{2};
                                L = varargin{3};
                                break
                            end
                        else
                            msg = 'Unable to parse argument ''%s''.';
                            error(msg,inputname(j+1))
                        end
                    end
                    
                otherwise
                    msg = 'Expected ''%s'' to be of class polygon or rectmap.';
                    error(msg,inputname(1))
                    
            end % switch
            
            % Retrieve options
            opt = scmapopt(opt);
            
            % Get data for the low-level functions
            w = vertex(poly);
            n = length(w);
            beta = angle(poly) - 1;
            
            if isempty(z)
                % Request corners
                if isempty(corner)
                    msg{1} = 'Select the images of the corners of the rectangle.';
                    msg{2} = 'Go in counterclockwise order and select a long rectangle edge first.';
                    corner = scselect(w,beta,4,'Select corners',msg);
                end
                
                % Apply SCFIX to enforce solver rules
                [w,beta,corner] = scfix('r',w,beta,corner);
                poly = polygon(w,beta+1);
                
                % Solve parameter problem (always necessary)
                [z,c,L,qdata] = rparam(w,beta,corner,opt.InitialGuess,opt);
                
            else
                % Prevertices, etc. given. Renumber to conform
                [w,beta,z] = rcorners(w,beta,z);
                nqpts = ceil(-log10(opt.Tolerance));
                qdata = scqdata(beta,nqpts);
                poly = polygon(w,beta+1);
            end
            
            map = map@scmap(poly,opt);
            
            map.prevertex = z;
            map.constant = c;
            map.stripL = L;
            map.qdata = qdata;
            
            % Now fill in apparent accuracy
            map.accuracy = accuracy(map);
            
        end
    end
end
