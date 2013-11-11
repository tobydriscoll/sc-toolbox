function M = mtimes(M,c)
%   Scale the image of a map by a complex constant.

% Usually this will be invoked by a child object, which needs to adjust
% its own constant as well.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: mtimes.m 33 1998-06-29 22:35:40Z tad $

% May need to swap arguments
if isa(M,'double') & isa(c,'scmap')
  tmp = M;
  M = c;
  c = tmp;
end

M.polygon = c*M.polygon;
