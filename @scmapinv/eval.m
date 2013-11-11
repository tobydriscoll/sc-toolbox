function zp = eval(Mi,wp)
%EVAL Evaluate inverse SC map.
%   EVAL(MI,WP), where MI is an SCMAPINV object and WP is a vector of points
%   in the polygon of the map, returns the inverse image of WP under the
%   map (the forward image under MI).
%   
%   See also SCMAPINV, SCMAPINV/SUBSREF.

%   Copyright 1998 by Toby Driscoll.
%   $Id: eval.m 7 1998-05-10 04:37:19Z tad $

zp = evalinv(Mi.themap,wp);
