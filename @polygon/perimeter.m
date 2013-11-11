function L = perimeter(p)
%PERIMETER Perimeter of a polygon.

%   $Id: perimeter.m 243 2003-03-03 16:18:13Z driscoll $

if isinf(p)
  L = Inf;
else
  w = vertex(p);
  L = sum( abs( diff( w([1:end 1]) ) ) );
end