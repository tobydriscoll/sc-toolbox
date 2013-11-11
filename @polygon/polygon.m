function p = polygon(x,y,alpha)
%POLYGON Contruct polygon object.
%   POLYGON(W) constructs a polygon object whose vertices are specified
%   by the complex vector W. Cusps and cracks are allowed.
%   
%   POLYGON(X,Y) specifies the vertices with two real vectors.
%   
%   POLYGON(W,ALPHA) or POLYGON(X,Y,ALPHA) manually specifies the interior
%   angles at the vertices, divided by pi.
%   
%   POLYGON accepts unbounded polygons (vertices at infinity). However,
%   you must supply ALPHA, and the vertices must be in counterclockwise
%   order about the interior.
%   
%   See also POLYGON/ANGLE, POLYGON/PLOT.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: polygon.m 71 1999-06-11 10:34:21Z tad $


superiorto('double');

w = [];
if nargin < 3
  alpha = [];
end

if nargin == 0
elseif isempty(x)
elseif isa(x,'polygon')
  p = x;
  return
elseif ~isreal(x) | nargin == 1 | (any(isinf(x)) & nargin==2)
  % Vertices passed as a complex vector
  w = x(:);
  % If first point is repeated at the end, delete the second copy
  % Thanks to Mark Embree for bug fix.
  if abs(w(end) - w(1)) < 3*eps
    w(end) = [];
  end
  if nargin > 1
    alpha = y;
  end
else
  % Vertices passed as two real vectors
  w = x(:) + i*y(:);
  % If first point is repeated at the end, delete the second copy
  if abs(w(end) - w(1)) < 3*eps
    w(end) = [];
  end  
end

% Create a polygon with current angles
p0.vertex = w(:);
p0.angle = alpha(:);
p0 = class(p0,'polygon');

n = length(w);

% Now compute angles if needed
if (n > 0) 
  [alpha,isccw,index] = angle(p0);
  if ~isccw
    p0.vertex = flipud(p0.vertex);
    alpha = flipud(alpha);
  end
  p0.angle = alpha;
  if abs(index) > 1
    warning('Polygon is multiple-sheeted.')
  end
end
  
p = p0;

