function fprime = dederiv(zp,z,beta,c)
%DEDERIV Derivative of the exterior map.
%   DEDERIV(ZP,Z,BETA,C) returns the derivative at the points of ZP of the
%   Schwarz-Christoffel exterior map defined by Z, BETA, and C.
%   
%   See also DEPARAM, DEMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: dederiv.m 7 1998-05-10 04:37:19Z tad $

% Support old syntax
if nargin < 4
  c = 1;
end

z = z(:);
beta = [beta(:);-2];
zprow = zp(:).';
fprime = zeros(size(zp));

npts = length(zp(:));
terms = 1 - zprow(ones(length(z),1),:)./z(:,ones(npts,1));
terms(length(z)+1,:) = zprow;
fprime(:) = c*exp(sum(log(terms).*beta(:,ones(npts,1))));
