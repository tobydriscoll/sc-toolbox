function M = plus(M,a)
%   Add a constant to the map (i.e., translate its image).

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: plus.m 36 1998-06-29 23:14:51Z tad $

% May need to swap arguments
if isa(M,'double') & isa(a,'crrectmap')
  tmp = M;
  M = a;
  a = tmp;
end

if length(a)==1 & isa(a,'double')
  M.diskmap = M.diskmap + a;
else
  error('Addition is not defined for these operands.')
end
