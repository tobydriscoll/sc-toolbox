function fprime = stderiv(zp,z,beta,c,j)
%STDERIV Derivative of the strip map.
%   STDERIV(ZP,Z,BETA,C) returns the derivative at the points of ZP of
%   the Schwarz-Christoffel strip map defined by Z, BETA, and C.
%   
%   See also STPARAM, STMAP.
    
%   Copyright 1998 by Toby Driscoll.
%   $Id: stderiv.m 298 2009-09-15 14:36:37Z driscoll $

%   If a fifth argument j is supplied, the terms corresponding to z(j)
%   are normalized by abs(zp-z(j)).  This is for Gauss-Jacobi
%   quadrature.

% Support old syntax
if nargin < 4
  c = 1;
end

log2 = 0.69314718055994531;
fprime = zeros(size(zp));
zprow = zp(:).';
npts = length(zprow);

% Strip out infinite prevertices
if length(z)==length(beta)
  ends = find(isinf(z));
  theta = diff(beta(ends));
  if z(ends(1)) < 0
    theta = -theta;
  end
  z(ends) = [];
  beta(ends) = [];
  % Adjust singularity index if given
  if nargin > 4
    j = j - (j > ends(1)) - (j > ends(2));
  end
else
  error('Vector of prevertices must include +/-Inf entries')
end
zcol = z(:);
bcol = beta(:);
n = length(z);

terms = -pi/2*(zprow(ones(n,1),:) - zcol(:,ones(npts,1)));
lower = (~imag(z));
terms(lower,:) = -terms(lower,:);
rt = real(terms);
big = abs(rt) > 40;
if any(any(~big))
  terms(~big) = log(-i*sinh(terms(~big)));
end
terms(big) = sign(rt(big)).*(terms(big)-i*pi/2) - log2;
if nargin > 4
  if j > 0
    terms(j,:) = terms(j,:)-log(abs(zprow-z(j)));
  end
end
fprime(:) = c*exp(pi/2*theta*zprow + sum(terms.*bcol(:,ones(npts,1))));
