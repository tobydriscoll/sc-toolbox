function wp = stmap(zp,w,beta,z,c,qdat)
%STMAP  Schwarz-Christoffel strip map.
%   STMAP(ZP,W,BETA,Z,C,QDAT) computes the values of the Schwarz-
%   Christoffel strip map at the points in vector ZP. The arguments W,
%   BETA, Z, C, and QDAT are as in STPARAM. STMAP returns a vector the
%   same size as ZP.
%
%   STMAP(ZP,W,BETA,Z,C,TOL) uses quadrature data intended to give an
%   answer accurate to within TOL.
%   
%   STMAP(ZP,W,BETA,Z,C) uses a tolerance of 1e-8.
%
%   See also STPARAM, STPLOT, STINVMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: stmap.m 298 2009-09-15 14:36:37Z driscoll $

if isempty(zp)
  wp = [];
  return
end

import sctool.*
n = length(w);
w = w(:);
beta = beta(:);
z = z(:);

% Quadrature data
if nargin < 6
  qdat = scqdata(beta,8);
elseif length(qdat)==1
  qdat = scqdata(beta,max(ceil(-log10(qdat)),8));
end
tol = 10^(-size(qdat,1));
  
shape = size(zp);
zp = zp(:);
zprow = zp.';
p = length(zp);
wp = zeros(p,1);

% For each point in zp, find the nearest prevertex.
[dist,sing] = min( abs(zprow(ones(n,1),:) - z(:,ones(1,p))) );
sing = sing(:);				% indices of prevertices

% Screen out images of prevertices
vertex = (dist(:) < tol);
wp(vertex) = w(sing(vertex));
leftend = (isinf(zp) & (real(zp) < 0));
wp(leftend) = w(z == -Inf);
rightend = (isinf(zp) & (real(zp) > 0));
wp(rightend) = w(z == Inf);
vertex = vertex | leftend | rightend;

% "Bad" points are closest to a prevertex of infinity.
atinf = find(isinf(w(:))); 		% infinite vertices 
bad = ismember(sing,atinf) & ~vertex;

if any(bad)
  % We can't begin integrations at a preimage of infinity. We will pick the
  % next-closest qualified prevertex.
  zf = z;
  zf(isinf(w)) = Inf;
  [tmp,s] = min( abs(zprow(ones(n,1),bad) - zf(:,ones(sum(bad),1))) );
  sing(bad) = s;
  
  % Because we no longer integrate from the closest prevertex, we must go in
  % stages to maintain accuracy.
  mid1 = real(z(s)) + i/2;
  mid2 = real(zp(bad)) + i/2;
else
  bad = zeros(p,1);		% all clear
end

% zs = the starting singularities
zs = z(sing);
% ws = map(zs)
ws = w(sing);

% Compute the map directly at "normal" points.
normal = ~bad & ~vertex;
if any(normal)
  I = stquad(zs(normal),zp(normal),sing(normal),z,beta,qdat);
  wp(normal) = ws(normal) + c*I;
end

% Compute map at "bad" points in stages.
if any(bad)
  I1 = stquad(zs(bad),mid1,sing(bad),z,beta,qdat);
  I2 = stquadh(mid1,mid2,zeros(sum(bad),1),z,beta,qdat);
  I3 = -stquad(zp(bad),mid2,zeros(sum(bad),1),z,beta,qdat);
  wp(bad) = ws(bad) + c*(I1 + I2 + I3);
end

wp = reshape(wp,shape);
