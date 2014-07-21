function beta = scangle(w)
%SCANGLE Turning angles of a polygon.
%   SCANGLE(W) computes the turning angles of the polygon whose vertices
%   are specified in the vector W.  The turning angle of a vertex
%   measures how much the heading changes at that vertex from the
%   incoming to the outgoing edge, normalized by pi.  For a finite
%   vertex, it is equal in absolute value to (exterior angle)/pi, with a
%   negative sign for left turns and positive for right turns.  Thus the
%   turn at a finite vertex is in (-1,1], with 1 meaning a slit.
%
%   At an infinite vertex the turning angle is in the range [-3,-1] and
%   is equal to the exterior angle of the two sides extended back from
%   infinity, minus 2.  SCANGLE cannot determine the angle at an
%   infinite vertex or its neighbors, and will return NaN's in those
%   positions.
%   
%   See also DRAWPOLY, DEMOINF.

%   Copyright 1998 by Toby Driscoll.
%   $Id: scangle.m 298 2009-09-15 14:36:37Z driscoll $

w = w(:);
n = length(w);
if n==0
  beta = [];
  return
end

atinf = isinf(w);
% These can't be determined.
mask = ~(atinf | atinf([2:n,1]) | atinf([n,1:n-1]));

dw = diff( w([n 1:n]) );
dwshift = dw([2:n,1]);
beta = NaN*ones(size(w));
beta(mask) = angle( dw(mask).*conj(dwshift(mask)) )/pi;

% It's ill-posed to tell a point (outward) from a slit (inward). Since
% the latter is much more common and important, we'll be generous in
% giving it the tie.
mod = abs(beta+1) < 1e-12;
beta(mod) = ones(size(beta(mod)));
