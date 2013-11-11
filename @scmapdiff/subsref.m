function wp = subsref(Md,S)
%SUBSREF Evaluate differentiated SC map by subscript notation.
%   MD(ZP), where MD is an SCMAPDIFF object and ZP is a vector of
%   points in canonical domain of the map, returns the derivative of the
%   map used to create MD at the points ZP.
%   
%   This just a synonym for EVAL(MD,ZP).
%
%   See also SCMAPDIFF, SCMAPDIFF/EVAL.

%   Copyright 1998 by Toby Driscoll.
%   $Id: subsref.m 7 1998-05-10 04:37:19Z tad $

if length(S) == 1 & strcmp(S.type,'()')
  wp = evaldiff(Md.themap,S.subs{1});
else
  error('Only syntax for SCMAPDIFF is a single parenthesized subscript.')
end

  