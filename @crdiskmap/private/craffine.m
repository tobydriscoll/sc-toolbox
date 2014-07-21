function [aff,wn] = craffine(w,beta,cr,Q,tol)
%CRAFFINE Affine transformations for crossratio formulation.
%   CRAFFINE(W,BETA,CR,Q) computes an (n-3)x2 array. The CR formulation
%   has n-3 implied "raw" conformal maps based on the S-C integral, each
%   of which is an affine transformation away from the true target
%   polygon. The output array describes these transformations. An extra
%   argument can be used to specify the desired tolerance.
%       
%   Because of the overlapping nature of the quadrilaterals in the
%   polygon's decomposition, only one side of W need be specified in
%   order to define these transformations uniquely. If W is an n-vector
%   of NaN's with at least two adjacent vertices filled in, [AFF,W] =
%   CRAFFINE(W,BETA,CR,Q) will also return the full polygon thus
%   defined. This is useful for computing maps to a rectangle, for
%   example, where the aspect ratio--hence the side lengths--of the
%   target W are not known in advance.
       
%   Copyright 1998 by Toby Driscoll.

import sctool.*
n = length(beta);
if nargin < 5
  tol = 1e-8;
  if nargin < 4
    Q = crqgraph(w);
  end
end

nqpts = max(4,ceil(-log10(tol)));
qdat = scqdata(beta,nqpts);

wn = NaN*ones(n,1);
aff = NaN*ones(n-3,2);
rawimage = ones(4,n-3);

% Deduce the side we will start with
sidenum = min(find(~isnan(w) & ~isnan(w([2:n,1]))));
if isempty(sidenum)
  error('You must specify at least one side of the target polygon.')
end
s = [sidenum,rem(sidenum,n)+1];
side = w(s);

% Start by embedding for the quadrilateral containing the boundary edge
% given by sidenum

% Which edge # is it?
edgenum = find(Q.edge(1,:)==s(1) & Q.edge(2,:)==s(2));
if isempty(edgenum)
  edgenum = find(Q.edge(2,:)==s(1) & Q.edge(1,:)==s(2));
end
% Which quadrilateral?
quadnum = min(find(any(Q.qledge==edgenum(ones(4,1),ones(n-3,1)))));
% Do the embedding & calculate affine constants
z = crembed(cr,Q,quadnum);
idx = Q.qlvert(:,quadnum);
v = -crquad(z(idx),idx,z,beta,qdat);
wn(idx) = v;
y = [wn(s),ones(2,1)] \ side(:);
aff(quadnum,:) = y(:).';
wn(idx) = wn(idx)*aff(quadnum,1) + aff(quadnum,2);

% Set up stack and begin 

% Keep track of "raw" quadrilateral images as they are found
rawimage(:,quadnum) = v;
stack = NaN*zeros(n-3,1);
newnbrs = Q.adjacent(:,quadnum) & isnan(aff(:,1)); % neighbors of quadnum...
m = sum(newnbrs);
stack(1:m) = find(newnbrs);		% ...go on the stack first
% Keep track of who put a diagonal on the stack
origin = NaN*zeros(n-3,1);
origin(1:m) = quadnum*ones(m,1);
stackptr = m;
while stackptr > 0
  q = stack(stackptr);
  oldq = origin(stackptr);
  stackptr = stackptr-1;
  
  % Embedding & raw image
  idx = Q.qlvert(:,q);
  z = crembed(cr,Q,q);
  newv = -crquad(z(idx),idx,z,beta,qdat);
  rawimage(:,q) = newv;

  oldv = rawimage(:,oldq);
  oldidx = Q.qlvert(:,oldq);
  
  % Find new transformation by composing with the old one
  [ref,oldref] = find(idx(:,ones(4,1))==oldidx(:,ones(4,1))');
  done = zeros(4,1);
  done(ref) = ones(3,1);
  y = [newv(ref),ones(3,1)] \ oldv(oldref);
  aff(q,1) = y(1)*aff(oldq,1);
  aff(q,2) = y(2)*aff(oldq,1) + aff(oldq,2);
  wn(idx(~done)) = newv(~done)*aff(q,1) + aff(q,2);

  newnbrs = Q.adjacent(:,q) & isnan(aff(:,1));
  m = sum(newnbrs);
  stack(stackptr+(1:m)) = find(newnbrs);
  origin(stackptr+(1:m)) = q*ones(m,1);
  stackptr = stackptr + m;
  
end
