function fprime = rsderiv(zp,z,beta,zb,c)
%RSDERIV Derivative of the Riemann surface map.

%   RSDERIV(ZP,Z,BETA,ZB,C) returns the derivative at the points of ZP of
%   the Schwarz-Christoffel disk map defined by Z, BETA, and C.

%   Copyright 2002 by Toby Driscoll.
%   $Id: rsderiv.m 298 2009-09-15 14:36:37Z driscoll $

% Support old syntax
if nargin < 4
  c = 1;
end

z = z(:);
beta = beta(:);
zprow = zp(:).';
fprime = zeros(size(zp));

npts = length(zp(:));
ZP = zprow(ones(length(beta),1),:);
Z = z(:,ones(npts,1));
terms = 1 - ZP./Z;
fp = c*exp(sum(log(terms).*beta(:,ones(npts,1))));

B = length(zb);
if B > 0
  ZB = zb(:,ones(1,npts));
  ZP = zprow(ones(B,1),:);
  fp = fp.*prod((ZP-ZB).*(1-ZP.*conj(ZB)),1);
end

fprime(:) = fp;