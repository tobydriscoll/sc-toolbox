function n = size(p,m)
%   Number of vertices.

%   Copyright 1998 by Toby Driscoll.
%   $Id: size.m 7 1998-05-10 04:37:19Z tad $

if nargin == 1
  n = [length(p.vertex) 1];
elseif m == 1
  n = length(p.vertex);
else
  n = 1;
end

  