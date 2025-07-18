function M = plus(M,a)
%   Add a constant to the map (i.e., translate its image).

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: plus.m 36 1998-06-29 23:14:51Z tad $

% May need to swap arguments
if isa(M,'double') & isa(a,'crdiskmap')
  tmp = M;
  M = a;
  a = tmp;
end

if isscalar(a) && isa(a,'double')
  M.affine(:,2) = M.affine(:,2) + a;
  M.scmap = M.scmap + a;
  M = center(M,M.center{1}+a);
else
  error('Addition is not defined for these operands.')
end
