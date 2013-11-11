function r = minus(p,q)
%   Translate a polygon, or subtract the vertices of two polygons.

%   Copyright 1999-2003 by Toby Driscoll.
%   $Id: minus.m 247 2003-03-03 16:28:22Z driscoll $

r = plus(p,-q);
