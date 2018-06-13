function [fprime,d2f] = dderiv(zp,z,beta,c)
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

%terms = 1 - zprow./z;
%fprime(:) = c*exp(sum(log(terms).*beta));
npts = length(zp(:));
terms = 1 - zprow(ones(length(beta),1),:)./z(:,ones(npts,1));
fprime(:) = c*exp(sum(log(terms).*beta(:,ones(npts,1))));

if nargout > 1   % 2nd derivative
    d2f = 0;
    for k = 1:length(z)
        d2f = d2f - (beta(k)/z(k))./terms(k,:).';
    end
    d2f = d2f .* fprime;
end

