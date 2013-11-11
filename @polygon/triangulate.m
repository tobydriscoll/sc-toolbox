function [tri,x,y] = triangulate(p,h)
%TRIANGULATE Triangulate the interior of a polygon.
%
%   [TRI,X,Y] = TRIANGULATE(P,H) attempts to find a reasonable
%   triangulation of the interior of polygon P so that the typical
%   triangle diameter is about H. If H is not specified, an automatic
%   choice is made.
%
%   If P is unbounded, the polygon is first truncated to fit in a
%   square.
%
%   TRIANGULATE uses DELAUNAY from Qhull, and as such does not have
%   guaranteed success for nonconvex regions. However, things should
%   go OK unless P has slits. 
%
%   See also TRUNCATE, DELAUNAY.

%   Copyright 2002 by Toby Driscoll.
%   $Id: triangulate.m 191 2002-09-10 18:21:26Z driscoll $
  
if isinf(p)
  warning('Truncating an unbounded polygon.')
  p = truncate(p);
end

w = vertex(p);
n = length(w);

if nargin < 2
  h = diam(p) / 40;
end

% Find points around boundary.
[wb,idx] = linspace(p,h/2);   % smaller spacing to help qhull

% On sides of a slit, put on extra points and perturb inward a little.
isslit = ( abs(angle(p)-2) < 10*eps );
slit = find( isslit | isslit([2:n 1]) );
[wbfine,idxfine] = linspace(p,h/6);
for k = slit(:)'
  new = (idxfine==k);
  old = (idx==k);
  wb = [ wb(~old); wbfine(new) ];  idx = [ idx(~old); idxfine(new) ];
  move = find(idx==k);
  normal = i*( w(rem(k,n)+1) - w(k) );
  wb(move) = wb(move) + 1e-8*normal;
end

% Take out points that are fairly close to a singularity, because it's
% hard to find the inverse mapping there.
for k = find(angle(p)<1)'
  close = abs(wb - w(k)) < h/3;
  wb(close) = [];  idx(close) = [];
end

% Add the polygon vertices.
wb = [ wb; w ];
idx = [ idx; (1:n)'];

% Find a hex pattern that covers the interior.
xlim = [ min(real(w)) max(real(w)) ];
ylim = [ min(imag(w)) max(imag(w)) ];
x = linspace(xlim(1),xlim(2),ceil(diff(xlim))/h+1);
y = linspace(ylim(1),ylim(2),ceil(diff(ylim))/h+1);
[X,Y] = meshgrid(x(2:end-1),y(2:end-1));
X(2:2:end,:) = X(2:2:end,:) + (x(2)-x(1))/2;

inside = isinpoly(X+1i*Y,p);
x = [ real(wb); X(inside) ];
y = [ imag(wb); Y(inside) ];

% Triangulate using qhull.
tri = delaunay(x,y);

% Those with a boundary vertex must be examined.
nb = length(wb);
check = find( any( tri<=nb, 2 ) );

% First, check triangle midpoints.
idx = tri(check,:);
z = x(idx) + 1i*y(idx);
out = ~isinpoly( sum(z,2)/3, p );

% On the rest, look for edges that cross two slit sides. 
check2 = find(~out);
sect1 = intersect(p,z(check2,[1 2]),1e-6);
sect2 = intersect(p,z(check2,[2 3]),1e-6);
sect3 = intersect(p,z(check2,[3 1]),1e-6);
out(check2( sum(sect1,2) > 1 )) = 1;
out(check2( sum(sect2,2) > 1 )) = 1;
out(check2( sum(sect3,2) > 1 )) = 1;

tri(check(out),:) = [];


