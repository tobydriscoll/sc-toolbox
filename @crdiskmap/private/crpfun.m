function f = crpfun(x,fdat)
%CRPFUN (not intended for calling directly by the user)
%   Nonlinear function for CRPARAM.
%   Copyright 1999 by Toby Driscoll.
%   $Id: crpfun.m 64 1999-01-29 00:58:06Z tad $

[n,beta,crtarget,Q,qdat] = deal(fdat{:});

crprever = exp(x);			% prevertex crossratios
crimage = zeros(n-3,1);			% image vertex crossratios

% Compute crossratio for each image quadrilateral
for k = 1:n-3
  prever = crembed(crprever,Q,k);
  w = -crquad(prever(Q.qlvert(:,k)),Q.qlvert(:,k),prever,beta,qdat);
  crimage(k) = (w(2)-w(1))*(w(4)-w(3))/((w(3)-w(2))*(w(1)-w(4)));
end

% Logarithmic scaling for residual
f = log(abs(crimage./crtarget));
