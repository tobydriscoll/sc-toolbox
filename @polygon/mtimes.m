function r = mtimes(p,q)
%   Multiplication of a polygon by a scalar.

%   Copyright 1998 by Toby Driscoll.
%   $Id: mtimes.m 7 1998-05-10 04:37:19Z tad $

if isa(q,'polygon')
  if isa(p,'polygon')
    error('Function ''*'' not defined for two polygon objects.')
  end
  tmp = p;
  p = q;
  q = tmp;
end

r = p;
r.vertex = r.vertex*q;
