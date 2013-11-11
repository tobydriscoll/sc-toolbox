function zp = subsref(Mi,S)
%SUBSREF Evaluate inverse map by subscript notation.
%   MI(WP), where MI is an SCMAPINV object and WP is a vector of points in
%   the polygon of the map, returns the inverse image of WP under the map.
%   
%   This just a synonym for EVAL(MI,WP), or EVALINV(INV(MI),WP).
%   
%   See also SCMAPINV, SCMAPINV/EVAL.

%   Copyright 1998 by Toby Driscoll.
%   $Id: subsref.m 7 1998-05-10 04:37:19Z tad $

if length(S) == 1 & strcmp(S.type,'()')
  zp = evalinv(Mi.themap,S.subs{1});
else
  error('Only syntax for SCMAPINV is a single parenthesized subscript.')
end

  