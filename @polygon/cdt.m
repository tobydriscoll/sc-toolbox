function T = cdt(p)
%CDT    Constrained Delaunay triangulation of polygon vertices.
%   T = CDT(P) returns a structure representing a constrained Delaunay
%   triangulation of the n polygon vertices. T has the fields:
%   
%      T.edge    : 2x(2n-3) matrix of indices of edge endpoints
%      T.triedge : 3x(n-2) matrix of triangle edge indices
%      T.edgetri : 2x(2n-3) matrix of triangle membership indices for
%                  the edges (boundary edges have a zero in 2nd row)
%   
%   See also PLOTCDT.

%   Copyright 1998 by Toby Driscoll.
%   $Id: cdt.m 7 1998-05-10 04:37:19Z tad $

w = p.vertex;
if any(isinf(w))
  error('CDT not possible for unbounded polygons.')
end

[e,te,et] = crtriang(w);
[e,te,et] = crcdt(w,e,te,et);

T = struct('edge',e,'triedge',te,'edgetri',et);
