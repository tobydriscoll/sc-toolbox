function M = plus(M,a)
%   Add a constant to the image of an SC map (i.e., translate image).

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: plus.m 32 1998-06-29 22:26:36Z tad $

% May need to swap arguments
if isa(M,'double') & isa(a,'scmap')
  tmp = M;
  M = a;
  a = tmp;
end

if length(a)==1 & isa(a,'double')
  M.polygon = M.polygon + a;
else
  error('Addition is not defined for these operands.')
end
