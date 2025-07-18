function fp = crderiv(zp,beta,cr,aff,wcfix,Q)
%CRDERIV  Derivative of the disk map in crossratio formulation.
%   CRDERIV(ZP,BETA,CR,AFF,WCFIX,Q) returns the derivative at the points
%   ZP of the Schwarz-Christoffel crossratio disk map. The arguments are
%   returned from CRPARAM, CRAFFINE, and CRFIXWC.
%
%   CRDERIV(ZP,BETA,CR,AFF,WCFIX,Q,TOL) uses quadrature data intended to
%   give an answer accurate to within roughly TOL.
%       
%   See also CRPARAM, CRAFFINE, CRFIXWC, CRMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crderiv.m 7 1998-05-10 04:37:19Z tad $

% Parse input and initialize
beta = beta(:);
fp = zeros(size(zp));
zp = zp(:).';

% Transform points into all embeddings, from the reference in wcfix
quadnum = wcfix(1);
mt = wcfix(2:5);
zl = (mt(1)*zp + mt(2)) ./ (mt(3)*zp + mt(4));
d0 = (mt(1)*mt(4) - mt(2)*mt(3)) ./ (mt(3)*zp + mt(4)).^2;
[zl,dl] = crspread(zl,quadnum,cr,Q);

% Choose best embeddings based on proximity to origin
[~,idx] = min(abs(zl));

% Compute derivatives via embeddings
for q = unique(idx)
  z = crembed(cr,Q,q);
  mask = (idx==q);
  % Compose Moebius transformations with disk map, and apply the affine
  % transformation constant
  fp(mask) = aff(q,1)*d0(mask).*dl(q,mask).*dderiv(zl(q,mask),z,beta);
end
