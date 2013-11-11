function M = plus(M,a)
%   Add a constant to the image of a map (i.e., translate image).

%   Copyright (c) 2002 by Toby Driscoll.
%   $Id: plus.m 206 2002-09-13 20:21:47Z driscoll $

% May need to swap arguments
if isa(M,'double') & isa(a,'riesurfmap')
  tmp = M;
  M = a;
  a = tmp;
end

if length(a)==1 & isa(a,'double')
  M.center = M.center + a;
  M.scmap = M.scmap + a;
else
  error('Addition is not defined for these operands.')
end
