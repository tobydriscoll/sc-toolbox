function M = mtimes(M,c)
%   Scale the image of a map by a complex constant.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: mtimes.m 33 1998-06-29 22:35:40Z tad $

% May need to swap arguments
if isa(M,'double') & isa(c,'stripmap')
  tmp = M;
  M = c;
  c = tmp;
end

M.constant = c*M.constant;
M.scmap = c*M.scmap;
