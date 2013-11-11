function M = mtimes(M,c)
%   Scale the image of a map by a complex constant.

%   Copyright 2002 by Toby Driscoll.
%   $Id: mtimes.m 204 2002-09-13 20:09:50Z driscoll $

% May need to swap arguments
if isa(M,'double') & isa(c,'riesurfmap')
  tmp = M;
  M = c;
  c = tmp;
end

M.constant = c*M.constant;
M.scmap = c*M.scmap;
M.center = c*M.center;