function [alpha,isccw,index] = angle(p)
%ANGLE Normalized interior angles of a polygon.
%   ALPHA = ANGLE(P) returns the interior angles, normalized by pi, of
%   the polygon P. 0 < ALPHA(J) <= 2 if vertex J is finite, and -2 <=
%   ALPHA(J) <= 0 if J is an infinite vertex. (This is consistent with
%   the definition as the angle swept from the exiting side through the
%   interior to the incoming side, with the vertices in counterclockwise
%   order.) It is impossible to compute angles for an unbounded polygon;
%   they must be supplied to the POLYGON constructor to be well-defined.
%   
%   See also POLYGON/POLYGON.

%   Copyright 1998 by Toby Driscoll.
%   $Id: angle.m 257 2003-03-27 15:48:46Z driscoll $

w = p.vertex; 
n = length(w);

if ~isempty(p.angle)
  % If angles have been assigned, return them
  alpha = p.angle;
else
  if isempty(w)
    alpha = [];
    isccw = [];
    index = [];
    return
  end
  
  if any(isinf(w))
    error('Cannot compute angles for unbounded polygons.')
  end

  % Compute angles 
  incoming = w - w([n 1:n-1]);
  outgoing = incoming([2:n,1]);
  alpha = mod( angle(-incoming.*conj(outgoing))/pi ,2 );

  % It's ill-posed to determine locally the slits (inward-pointing) from
  % points (outward-pointing). Check suspicious cases at the tips to see if
  % they are interior to the rest of the polygon.
  mask = (alpha < 100*eps) | (2-alpha < 100*eps);
  if all(mask)
    % This can happen if all vertices are collinear
    alpha(:) = 0;
    isccw = 1;				% irrelevant
    index = 1;                          % irrelevant
    return
  end
  slit = logical(isinpoly(w(mask),w(~mask)));
  fmask = find(mask);
  alpha(fmask(slit)) = 2;
  alpha(fmask(~slit)) = 0;
  
end

% Now test--if incorrect, assume the orientation is clockwise
index = sum(alpha-1)/2;                 % should be integer
if abs(index - round(index)) > 100*sqrt(n)*eps
  % Try reversing the interpretation of a crack
  mask = (alpha < 2*eps) | (2-alpha < 2*eps);
  alpha(~mask) = 2 - alpha(~mask);
  index = sum(alpha-1)/2;                 % should be integer
  % If still not OK, something is wrong
  if abs(index - round(index) ) > 100*sqrt(n)*eps
    error('Invalid polygon.')
  end
end

index = round(index);
isccw = (index < 0);
    
