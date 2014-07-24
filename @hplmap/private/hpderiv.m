function fprime = hpderiv(zp,z,beta,c)
%HPDERIV Derivative of the half-plane map.
%   HPDERIV(ZP,Z,BETA,C) returns the derivative at the points of ZP of
%   the Schwarz-Christoffel half-plane map defined by Z, BETA, and C.
%
%   See also HPPARAM, HPMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: hpderiv.m 298 2009-09-15 14:36:37Z driscoll $

% Support old syntax
if nargin < 4
  c = 1;
end

zf = z(~isinf(z));
beta = beta(~isinf(z));
zprow = zp(:).';
fprime = zeros(size(zp));

npts = length(zp(:));
terms = zprow(ones(length(beta),1),:) - zf(:,ones(npts,1));
fprime(:) = c*exp(sum(log(terms).*beta(:,ones(npts,1))));
