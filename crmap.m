function wp = crmap(zp,w,beta,cr,aff,wcfix,Q,qdat)
%CRMAP  Schwarz-Christoffel disk map in crossratio formulation.
%   CRMAP(ZP,W,BETA,CR,AFF,WCFIX,Q,QDAT) computes the values of the disk
%   map at the points in vector ZP. The arguments are returned from
%   CRPARAM, CRAFFINE, and CRFIXWC.
%
%   CRMAP(ZP,W,BETA,CR,AFF,WCFIX,Q,TOL) uses quadrature data intended to
%   give an answer accurate to within roughly TOL.
%       
%   CRMAP(ZP,W,BETA,CR,AFF,WCFIX,Q) uses a tolerance of 1e-8.
%
%   See also CRPARAM, CRAFFINE, CRFIXWC, CRPLOT, CRINVMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crmap.m 7 1998-05-10 04:37:19Z tad $

% Parse input and initialize
n = length(w);
w = w(:);
beta = beta(:);
if nargin < 8
  qdat = scqdata(beta,8);
elseif length(qdat)==1
  qdat = scqdata(beta,max(ceil(-log10(qdat)),4));
end
wp = zeros(size(zp));
zp = zp(:);
p = length(zp);

% Transform points into all embeddings, from the reference in wcfix
quadnum = wcfix(1);
zl = (wcfix(2)*zp + wcfix(3))./(wcfix(4)*zp + wcfix(5));
zl = crspread(zl,quadnum,cr,Q);

% Choose best embeddings based on proximity to origin
[tmp,idx] = min(abs(zl));

% Compute maps via embeddings
for q = unique(idx)
  z = crembed(cr,Q,q);
  mask = (idx==q);
  wp(mask) = crmap0(zl(q,mask),z,beta,aff(q,:),qdat);
end
