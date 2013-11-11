function idx = isinpoly(wp,p,varargin)
%ISINPOLY Identify points interior/exterior to a polygon.
%   ISINPOLY(WP,P) returns a logical vector the size of WP in which
%   nonzero means the corresponding point is inside polygon P and zero
%   means it is outside. 
%
%   ISINPOLY(WP,P,TOL) considers points within TOL of the boundary to be
%   inside P. Without this argument, points on the boundary may or may not
%   register as inside.
%
%   See also POLYGON/WINDING.

%   Copyright 1998-2003 by Toby Driscoll.
%   $Id: isinpoly.m 231 2003-01-09 14:49:49Z driscoll $

idx = logical( winding(p,wp,varargin{:}) );

