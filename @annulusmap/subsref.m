function wp = subsref(M,S)
%SUBSREF Evaluate map by subscript notation.
%   M(ZP), where M is a SC map and ZP is a vector of points in the
%   canonical domain of the map, returns the image of the points in ZP. 
%   
%   This just a synonym for EVAL(M,ZP).
%   
%   See also EVAL, SCMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: subsref.m 7 1998-05-10 04:37:19Z tad $

if length(S) == 1 & strcmp(S.type,'()')
  wp = eval(M,S.subs{1});
else
  error('Only syntax for SCMAP is a single parenthesized subscript.')
end

  