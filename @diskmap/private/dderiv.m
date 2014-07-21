function fprime = dderiv(zp,z,beta,c)
%DDERIV Derivative of the disk map.
%   DDERIV(ZP,Z,BETA,C) returns the derivative at the points of ZP of
%   the Schwarz-Christoffel disk map defined by Z, BETA, and C.
%
%   See also DPARAM, DMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: dderiv.m 7 1998-05-10 04:37:19Z tad $

% Support old syntax
if nargin < 4
  c = 1;
end

z = z(:);
beta = beta(:);
zprow = zp(:).';
fprime = zeros(size(zp));

npts = length(zp(:));
terms = 1 - zprow(ones(length(beta),1),:)./z(:,ones(npts,1));
fprime(:) = c*exp(sum(log(terms).*beta(:,ones(npts,1))));
