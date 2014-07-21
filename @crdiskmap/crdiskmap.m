classdef (InferiorClasses = {?double}) crdiskmap < scmap
    
    properties
        crossratio
        affine
        qlgraph
        qdata
        original
        center
        accuracy
    end
    
    methods
        function map = crdiskmap(poly,varargin)
            %CRDISKMAP Schwarz-Christoffel cross-ratio disk map object.
            %   CRDISKMAP(P) constructs a Schwarz-Christoffel crossratio disk map
            %   object for the polygon P. The parameter problem is solved using
            %   default options for the crossratios of the prevertices.
            %
            %   CRDISKMAP(P,OPTIONS) uses an options structure of the type created
            %   by SCMAPOPT in solving the parameter problem.
            %
            %   CRDISKMAP(P,CR,Q) creates a crdiskmap object having the given
            %   prevertex crossratios CR and quadrilateral graph Q.  An OPTIONS
            %   argument can be added, although only the error tolerance will be
            %   used.
            %
            %   CRDISKMAP(M), where M is a crdiskmap object, just returns M.
            %
            %   CRDISKMAP(M,P) returns a new crdiskmap object for the polygon P
            %   using the options in crdiskmap M. The crossratios in M will be used
            %   as the starting guess for the parameter problem of the new map. Thus
            %   P should properly be a perturbation of the polygon for M. An OPTIONS
            %   structure may be added to override options in M.
            %
            %   Use CRDISKMAP instead of DISKMAP when the polygon has elongations,
            %   or when DISKMAP fails to converge.
            %
            %   See also SCMAPOPT, classes POLYGON, DISKMAP.
            
            %   Copyright 1998-2001 by Toby Driscoll.
            
            % Assign empties to optional args
            import sctool.*
            cr = [];
            c = [];
            opt = [];
            
            % Branch based on class of first argument
            switch class(poly)
                
                case 'crdiskmap'
                    oldmap = poly;
                    % Continuation of given map to given polygon
                    poly = varargin{1};
                    opt = scmapopt(oldmap);
                    cr0 = oldmap.crossratio;
                    if length(cr0) ~= length(poly)-3
                        msg = 'Polygon %s must have the same length as that in %s.';
                        error(msg,inputname(2),inputname(1))
                    end
                    if nargin > 2
                        opt = scmapopt(opt,varargin{2});
                    end
                    opt = scmapopt(opt,'initial',cr0);
                    orig = oldmap.original;
                    
                case 'polygon'
                    % Parse optional arguments
                    for j = 1:length(varargin)
                        arg = varargin{j};
                        % Each arg is an options struct, cr, or quadrilateral graph
                        if isa(arg,'struct')
                            if strmatch('edge',fieldnames(arg))
                                Q = arg;
                            else
                                opt = arg;
                            end
                        elseif length(arg) == length(poly)-3
                            cr = arg;
                            cr = cr(:);
                        else
                            msg = 'Unable to parse argument ''%s''.';
                            error(msg,inputname(j+1))
                        end
                    end
                    
                otherwise
                    msg = 'Expected ''%s'' to be of class polygon or crdiskmap.';
                    error(msg,inputname(1))
                    
            end % switch
            
            
            % Retrieve options
            opt = scmapopt(opt);
            
            % Get data for the low-level functions
            w = vertex(poly);
            beta = angle(poly)-1;
            
            % Find prevertices if necessary
            if isempty(cr)
                % Apply SCFIX to enforce solver rules
                [w,beta] = scfix('d',w,beta);
                
                % Solve
                if isempty(opt.InitialGuess)
                    % Standard solution
                    [w,beta,cr,aff,Q,orig,qdata] = ...
                        crparam(w,beta,opt.InitialGuess,opt);
                    % Remake polygon to reflect change in w
                    poly = polygon(w,beta+1);
                else
                    % Continutation syntax
                    [cr,aff,Q] = crparam(w,beta,opt.InitialGuess,opt);
                    nqpts = ceil(-log10(opt.Tolerance));
                    qdata = scqdata(beta,nqpts);
                end
                
            else
                orig = true(size(w));
                % Base accuracy of quadrature on given options
                nqpts = ceil(-log10(opt.Tolerance));
                qdata = scqdata(beta,nqpts);
            end
            
            map = map@scmap(poly,opt);
            
            map.crossratio = cr;
            map.affine = aff;
            map.qlgraph = Q;
            map.qdata = qdata;
            map.original = orig;
            
            % Set conformal center as center of 1st triangle
            T = Q.qlvert(1:3,1);
            wc = mean(w(T));
            wcfix = crfixwc(w,beta,cr,aff,Q,wc);
            map.center = {wc,wcfix};
                       
            % Now fill in true accuracy
            map.accuracy = accuracy(map);
            
        end
    end
    
    methods (Static)
        % Needed by crrectmap
        function varargout = craffine(varargin)
            [varargout{1:nargout}] = craffine(varargin{:});
        end
        function varargout = crembed(varargin)
            [varargout{1:nargout}] = crembed(varargin{:});
        end
        function varargout = imap0(varargin)
            [varargout{1:nargout}] = crimap0(varargin{:});
        end
        function varargout = map0(varargin)
            [varargout{1:nargout}] = crmap0(varargin{:});
        end

    end
end
