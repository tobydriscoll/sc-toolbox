function [u,uquad,zr] = crgather(u,uquad,quadnum,cr,Q,zr)
%CRGATHER Convert points into a single embedding in CR formulation.
%   Each quadrilateral (crossratio) has an associated embedding of the
%   prevertices. These embeddings are linked by Moebius transformations,
%   each of which is well-conditioned.  CRGATHER(U,UQUAD,QUADNUM,CR,Q)
%   assumes that the points of U are given in the embeddings described
%   by UQUAD. The Moebius transformations are applied recursively to map
%   the points to their representation in the single embedding QUADNUM.
%       
%   See also CRPARAM, CRSPREAD, MOEBIUS.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crgather.m 42 1998-07-01 19:02:00Z tad $

n3 = length(cr);

if nargin < 6
  % Initial call (nonrecursive)
  zr = NaN*zeros(4,n3);
end

% Place the quadrilateral prevertices in a rectangle around the origin.
idx = Q.qlvert(:,quadnum);
r = cr(quadnum);
f1 = sqrt(1/(r+1));
f2 = sqrt(r/(r+1));
zr(:,quadnum) = f1*[-1;-1;1;1] + i*f2*[-1;1;1;-1];

% Recurse on neighbors to map into embedding quadnum
nbr = Q.adjacent(quadnum,:) & isnan(zr(1,:));
for q = find(nbr)
  % First, map into embedding q
  [u,uquad,zr] = crgather(u,uquad,q,cr,Q,zr);
  % Find the 3 points in common with q
  [i1,i2] = find(Q.qlvert(:,quadnum*ones(4,1))==Q.qlvert(:,q*ones(4,1))');
  % Map from q to quadnum
  mt = double(moebius(zr(i2,q),zr(i1,quadnum)));
  mask = (uquad==q);
  if any(mask)
    u(mask) = (mt(2)*u(mask) + mt(1))./(mt(4)*u(mask) + mt(3));
    uquad(mask) = quadnum*ones(sum(mask),1);
  end
end
