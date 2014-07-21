function fp = crrderiv(zp,w,beta,wr,betar,cr,aff,affr,Q,qdat,qdatr)
%CRRDERIV  Derivative of the crossratio rectified map.
%   CRRDERIV(ZP,W,BETA,WR,BETAR,CR,AFF,AFFR,Q) returns the derivative at
%   the points ZP of the Schwarz-Christoffel crossratio rectified
%   map. The arguments are returned from CRPARAM, CRAFFINE, and CRRECT.
%
%   CRRDERIV(ZP,W,BETA,WR,BETAR,CR,AFF,AFFR,Q,TOL) uses quadrature data
%   intended to give an answer accurate to within roughly TOL.
%       
%   See also CRPARAM, CRAFFINE, CRRECT, CRRMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crrderiv.m 89 1999-11-23 18:32:52Z tad $

% Parse input and initialize
beta = beta(:);
n = length(beta);
if nargin < 11
  if nargin < 10 
    qd = cell(0,0);
  else
    qd{1} = qdat;
  end
else
  qd = {qdat,qdatr};
end

fp = zeros(size(zp));

% Compute inverse images in disk
[wp,qnum,up] = crrmap(zp,w,beta,wr,betar,cr,aff,affr,Q,qd{:});
qnum = qnum(:).';

% Compute derivatives via embeddings
for q = unique(qnum(~isnan(qnum)))
  mask = (qnum==q);
  z = crembed(cr,Q,q);
  % Compose disk map with the affine transformation constant
  fp(mask) = aff(q,1)*dderiv(up(mask),z,beta);
  % Do same for inverse of rectified
  fp(mask) = fp(mask)./(affr(q,1)*dderiv(up(mask),z,betar));
end
