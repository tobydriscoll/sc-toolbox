function [index,onvtx] = isinpoly(z,w,beta,tol)
%ISINPOLY Identify points inside a polygon.
%   ISINPOLY(z,w) or ISINPOLY(z,w,beta) returns a vector the size of z
%   such that each nonzero corresponds to a point inside the polygon
%   defined by w and beta.
%   
%   More precisely, the value returned for a point is the winding number
%   of the polygon about that point.
%   
%   The problem becomes ill-defined for points very near an edge or
%   vertex. ISINPOLY(z,w,tol) or ISINPOLY(z,w,beta,tol) considers points
%   within roughly tol of the boundary to be "inside", and computes
%   winding number for such points as the number of conformal images the
%   point ought to have.

%   Copyright 1998 by Toby Driscoll.
%   $Id: isinpoly.m 298 2009-09-15 14:36:37Z driscoll $

% Uses the argument principle, with some gymnastics for boundary points. 
import sctool.*

if nargin < 4
  tol = eps;
  if nargin < 3
    beta = scangle(w);
  else
    if length(beta)==1
      tol = beta;
      beta = scangle(w);
    end
  end
end

index = zeros(size(z));
n = length(w);

% Rescale to make differences relative
scale = mean(abs(diff(w([1:n 1]))));

% Trivial case (e.g., a single or repeated point)
if ~any(scale > eps)
  return
end

w = w/scale;
z = z/scale;

% Array of differences between each z and each w
zr = z(:).';
np = length(zr);
d = w(:,ones(np,1)) - zr(ones(n,1),:);

% Avoid divides by zero
bad = abs(d) < eps;
d(bad) = eps*ones(sum(bad(:)),1);

% Diffs of imag(log(w-z)) around the polygon 
ang = angle(d([2:n,1],:)./d) / pi;

% Find boundary points (edge and vertex)
tangents = sign(w([2:n 1]) - w);

% If points are repeated (e.g. crowding), skip to new point
for p = find( tangents.' == 0 );
  v = [w(p+1:n);w(1:n)];
  g = find(v ~= w(p));
  tangents(p) = sign(v(g(1)) - w(p));
end
  
% Points which are close to an edge
onbdy = (abs(imag(d./tangents(:,ones(np,1)))) < 10*tol);
% Points which are essentially vertices
onvtx = (abs(d) < tol);
% Correction: points must be on the closed, finite edge segment
onbdy = onbdy & ( (abs(ang) > .9) | onvtx | onvtx([2:n 1],:) );

% Truly interior points are easy: add up the args
interior = ~any(onbdy);
if any(interior)
  index(interior) = round(sum(ang(:,interior))/2);
end

% Boundary points are tricky
for k = find(~interior)
  % Index wrt other parts of polygon
  S = sum(ang(~onbdy(:,k),k));
  % We pretend a vertex point is on either adjacent side
  b = beta(onvtx(:,k));
  % Each edge membership counts as 1/2 winding number (either sign)
  augment = sum(onbdy(:,k)) - sum(onvtx(:,k)) - sum(b);

  index(k) = round(augment*sign(S) + S)/2;
end

index = logical(index);