function zp = crinvmap(wp,w,beta,cr,aff,wcfix,Q,qdat,options)
%CRINVMAP S-C disk inverse map in crossratio formulation.
%   CRINVMAP(WP,W,BETA,CR,AFF,WCFIX,Q) computes the inverse of the disk
%   map with given conformal center. You may append the optional
%   parameters QDAT and OPTIONS as in DINVMAP.
%       
%   You must first run CRPARAM, CRAFFINE, and CRFIXWC.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crinvmap.m 88 1999-11-23 18:20:45Z tad $

% Parse input and initialize
import sctool.*
n = length(w);
beta = beta(:);
zp = zeros(size(wp));
wp = wp(:);
lenwp = length(wp);
if nargin < 9
  options = [];
  if nargin < 8
    qdat = [];
  end
end

if isempty(qdat)
  qdat = scqdata(beta,8);
elseif length(qdat)==1
  qdat = scqdata(beta,max(ceil(-log10(qdat)),2));
end

% For each embedding, perform inverse maps for appropriate points
quadnum = zeros(lenwp,1);		% keep track of embeddings
for q=1:n-3
  idx = find(~quadnum);
  mask = abs(isinpoly(wp(idx),w(Q.qlvert(:,q)),10^(-size(qdat,1))));
  if any(mask)
    idx = idx(logical(mask));
    z = crembed(cr,Q,q);
    zp(idx) = crimap0(wp(idx),z,beta,aff(q,:),qdat,options);
    quadnum(idx) = q*ones(length(idx),1);
  end
  if all(quadnum), break, end
end

% Convert from local embeddings to global one
zp = crgather(zp,quadnum,wcfix(1),cr,Q);
mt = wcfix(2:5);
zp = (-mt(4)*zp + mt(2))./(mt(3)*zp - mt(1));
