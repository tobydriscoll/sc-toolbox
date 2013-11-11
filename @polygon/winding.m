function idx = winding(p,wp,varargin)
%WINDING Winding number of points with respect to a polygon.
%   WINDING(P,WP) returns a vector the size of WP of winding numbers with
%   respect to P. A zero value means the point is outside P; a value
%   greater than 1 means it lies on multiple sheets. 
%
%   WINDING(P,WP,TOL) makes the boundary of P "fuzzy" by a distance
%   TOL. This may be needed to compute winding number for points on the
%   boundary that you want to be considered "just inside." 
%
%   See also POLYGON/ISINPOLY.

%   Copyright 2003 by Toby Driscoll.
%   $Id: winding.m 230 2003-01-09 14:48:25Z driscoll $

if isinf(p)
  warning('SC:Truncation','Using a truncated version of the polygon.')
  p = truncate(p);
end

idx = double( isinpoly(wp,p.vertex,varargin{:}) );
