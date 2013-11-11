function M = crdiskmap(poly,varargin)
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
%   $Id: crdiskmap.m 129 2001-05-07 15:04:13Z driscoll $

superiorto('double');


if nargin == 0
  M.crossratio = [];
  M.affine = [];
  M.qlgraph = [];
  M.qdata = [];
  M.original = [];
  M.center = {};
  M_sc = scmap;
  M.accuracy = [];
  M = class(M,'crdiskmap',M_sc);
  return
end

% Assign empties to optional args
cr = [];
c = [];
opt = [];

% Branch based on class of first argument
switch class(poly)

case 'crdiskmap'
  M = poly;
  if nargin == 1
    % Self-return
    return
  else
    % Continuation of given map to given polygon
    poly = varargin{1};
    opt = scmapopt(M);
    cr0 = M.crossratio;
    if length(cr0) ~= length(poly)-3
      msg = 'Polygon %s must have the same length as that in %s.';
      error(sprintf(msg,inputname(2),inputname(1)))
    end
    if nargin > 2
      opt = scmapopt(opt,varargin{2});
    end
    opt = scmapopt(opt,'initial',cr0);
    orig = M.original;
  end
  
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
      error(sprintf(msg,inputname(j+1)))
    end
  end
  
otherwise
  msg = 'Expected ''%s'' to be of class polygon or crdiskmap.';
  error(sprintf(msg,inputname(1)))
  
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
  orig = logical(ones(size(w)));
  % Base accuracy of quadrature on given options
  nqpts = ceil(-log10(opt.Tolerance));
  qdata = scqdata(beta,nqpts);
end

M.crossratio = cr;
M.affine = aff;
M.qlgraph = Q;
M.qdata = qdata;
M.original = orig;

% Set conformal center as center of 1st triangle
T = Q.qlvert(1:3,1);
wc = mean(w(T));
wcfix = crfixwc(w,beta,cr,aff,Q,wc);
M.center = {wc,wcfix};

% Make a parent scmap object
M_sc = scmap(poly,opt);

% Leave a spot for accuracy and create object
M.accuracy = [];
if ~isa(M,'crdiskmap')
  M = class(M,'crdiskmap',M_sc);
else
  M.scmap = M_sc;
end

% Now fill in true accuracy
M.accuracy = accuracy(M);

