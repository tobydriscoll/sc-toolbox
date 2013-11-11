function M = crrectmap(poly,varargin)
%CRRECTMAP Schwarz-Christoffel cross-ratio disk map object.
%   CRRECTMAP(P) constructs a Schwarz-Christoffel crossratio rectified
%   map object for the polygon P. The parameter problem is solved, using
%   default options, for the crossratios of the prevertices on the
%   disk. Then you are required to graphically assign angles of the
%   rectified polygon. These angles determine the rectilinear polygon
%   that is considered the canonical domain.
%   
%   CRRECTMAP(P,OPTIONS) uses an options structure of the type created
%   by SCMAPOPT in solving the parameter problem.
%   
%   CRRECTMAP(P,ALPHAR) uses supplied angles for the rectified
%   polygon. Here ALPHAR is a vector having the same length as P. Each
%   element is an interior angle of the rectified polygon, normalized by
%   pi. The only valid entries are {.5,1,1.5,2}. An OPTIONS argument can
%   be added.
% 
%   CRRRECTMAP(M,ALPHAR) modifies the previously computed crrectmap M to
%   have rectified angles ALPHAR. Make ALPHAR the empty matrix to choose
%   new angles graphically. 
%   
%   CRRECTMAP(M,P), CRRECTMAP(M,P,ALPHAR) and CRRECTMAP(M,P,OPTIONS) use
%   continutation to go from the polygon in crrectmap M to P.
%   
%   CRRECTMAP(M), where M is a crrectmap object, just returns M.
%   
%   CRRECTMAP(MD), where MD is a crdiskmap object, uses the polygon and
%   prevertex crossratios already specified in MD, and requests angle
%   selection. CRRECTMAP(MD,ALPHAR) is also valid. 
%   
%   Use CRRECTMAP instead of RECTMAP when the polygon has multiple
%   elongations, or when RECTMAP fails to converge.
%   
%   See also SCMAPOPT, classes POLYGON, CRDISKMAP, RECTMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crrectmap.m 7 1998-05-10 04:37:19Z tad $

superiorto('double');

if nargin==0
  M.diskmap = crdiskmap;
  M.rectpolygon = [];
  M.rectaffine = [];
  M.rectqdata = [];
  M_sc = scmap;
  M = class(M,'crrectmap',M_sc);
  return
end

% Assign empties to optional args
alphar = [];
opt = [];
MD = [];

% Branch based on class of first argument
switch class(poly)

  case 'crrectmap'
    if nargin == 1
      % Self-return
      M = poly;
      return
    else
      M0 = poly;
      if isa(varargin{1},'double')
	% New alphar
	MD = M0.diskmap;
	alphar = varargin{1};
      elseif isa(varargin{1},'polygon')
	% Continuation
	MD = crdiskmap(M0.diskmap,varargin{1});
	if nargin > 2
	  if isstruct(varargin{2})
	    opt = varargin{2};
	    alphar = [];
	  else
	    alphar = varargin{2};
	  end
	else
	  alphar = angle(M0.rectpolygon);
	end
      else
	error('Invalid continuation syntax.')
      end
      poly = polygon(MD);
      opt = scmapopt(MD);
    end
  
  case 'crdiskmap'
    MD = poly;
    poly = polygon(MD);
    if nargin > 1
      alphar = varargin{1};
    end
    
  case 'polygon'
    % Parse optional arguments
    for j = 1:length(varargin)
      arg = varargin{j};
      % Each arg is an options struct or alphar
      if isa(arg,'struct')
	opt = arg;
      elseif length(arg) == length(poly)
	alphar = arg(:);
      else
	msg = 'Unable to parse argument ''%s''.';
      error(sprintf(msg,inputname(j+1)))
    end
  end
  
otherwise
  msg = 'Expected ''%s'' to be of class polygon, crrectmap, or crdiskmap.';
  error(sprintf(msg,inputname(1)))
  
end % switch


% Retrieve options
opt = scmapopt(opt);

% Get data for the low-level functions
w = vertex(poly);
beta = angle(poly) - 1;

% Find prevertex crossratios if necessary
if isempty(MD)
  MD = crdiskmap(poly,opt);
  param = parameters(MD);
  orig = param.original;
  poly = polygon(MD);
  w = vertex(poly);
  beta = angle(poly) - 1;
else
  param = parameters(MD);
  orig = logical(ones(size(w)));
  opt = scmapopt(MD);
end

M.diskmap = MD;

cr = param.crossratio;
aff = param.affine;
Q = param.qlgraph;

% Open rectified angle selection if necessary
if isempty(alphar)
  [wr,betar,affr] = crrect(w,beta,cr,aff,Q);
else
  % Check entries for validity
  if any(~ismember(alphar,[.5 1 1.5 2])) | abs(sum(alphar-1)+2) > 1e-12
    error('Invalid rectified angles supplied.')
  end
  % Convert alphas to betas, expanding if polygon was split
  betar = zeros(size(w));
  betar(orig) = alphar - 1;
  [wr,betar,affr] = crrect(w,beta,cr,aff,Q,betar);
end

polyr = polygon(wr,betar+1);
M.rectpolygon = polyr;
M.rectaffine = affr;
M.rectqdata = scqdata(betar,ceil(-log10(opt.Tolerance)));

% Make a parent scmap object
M_sc = scmap(poly,opt);

% Create object. Even though MD is a field, we're not really a child.
M = class(M,'crrectmap',M_sc);
