function fprime = rderiv(zp,z,beta,c,L,zs)
%RDERIV Derivative of the rectangle map.
%   RDERIV(ZP,Z,BETA,C,L) returns the derivative at the points of ZP of
%   the Schwarz-Christoffel rectangle map defined by Z, BETA, C, and L.
%   
%   If a sixth argument is supplied, it is assumed to be the image of Z
%   on the intermediate strip; see R2STRIP.
%   
%   See also RPARAM, RMAP, R2STRIP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: rderiv.m 298 2009-09-15 14:36:37Z driscoll $

n = length(z);

if nargin < 6
  % Find prevertices on the strip
  zs = r2strip(z,z,L);
  zs = real(zs) + i*round(imag(zs)); 	% put them *exactly* on edges
end

% First compute map and derivative from rectangle to strip
[F,dF] = r2strip(zp,z,L);

% Now compute derivative of map from strip to polygon
% Add in ends of strip
ends = find(diff(imag(z([1:n 1]))));
zs = [zs(1:ends(1));Inf;zs(ends(1)+1:ends(2));-Inf;zs(ends(2)+1:n)];
bs = [beta(1:ends(1));0;beta(ends(1)+1:ends(2));0;beta(ends(2)+1:n)];
dG = stripmap.deriv(F,zs,bs);

% Put it together
fprime = c*dF.*dG;

