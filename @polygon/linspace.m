function [z,idx] = linspace(p,m)
%LINSPACE Evenly spaced points around the polygon.
%   LINSPACE(P,N) returns a vector of N points evenly spaced on the
%   polygon P, starting with the first vertex.
%
%   LINSPACE(P,H) for H<1 instead uses H as an upper bound on the arc
%   length between points.
%
%   [Z,IDX] = LINSPACE(...) returns the points and an identically sized
%   vector of corresponding side indices.
%
%   If the polygon is unbounded, an error results.

%   Copyright 1998-2002 by Toby Driscoll. 
%   $Id: linspace.m 182 2002-09-05 15:40:11Z driscoll $

n = length(p);
w = vertex(p);
dw = diff(w([1:n 1]));

if any(isinf(w))
  error('LINSPACE cannot be applied to unbounded polygons.')
end

% Arc lengths of sides.
s = abs(dw);
s = cumsum([0;s]);
L = s(end);
s = s/L;  % relative arc length

% Evenly spaced points in arc length.
if m < 1
  % How many points will we need?
  m = ceil(L/m) + 1;
end
zs = (0:m-1)'/m;
z = zs;
done = logical(zeros(size(z)));
idx = zeros(size(z));

% Translate to polygon sides.
for j = 1:n
  mask = (~done) & (zs < s(j+1));
  z(mask) = w(j) + dw(j)*(zs(mask)-s(j))/(s(j+1)-s(j));
  idx(mask) = j;
  done = mask | done;
end

