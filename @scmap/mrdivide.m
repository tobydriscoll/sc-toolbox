function M = mrdivide(M,a)
%   Divide the image of an SC map by a constant.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: mrdivide.m 44 1998-07-01 20:21:54Z tad $

if isa(a,'double')
  M = M * (1/a);
else
  error('Cannot divide by an SC map.')
end

