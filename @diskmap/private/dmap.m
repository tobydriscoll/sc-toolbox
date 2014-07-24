function wp = dmap(zp,w,beta,z,c,qdat)
%DMAP  Schwarz-Christoffel disk map.
%   DMAP(ZP,W,BETA,Z,C,QDAT) computes the values of the Schwarz-
%   Christoffel disk map at the points in vector ZP. The arguments W,
%   BETA, Z, C, and QDAT are as in DPARAM. DMAP returns a vector the
%   same size as ZP.
%
%   DMAP(ZP,W,BETA,Z,C,TOL) uses quadrature data intended to give an
%   answer accurate to within TOL.
%   
%   DMAP(ZP,W,BETA,Z,C) uses a tolerance of 1e-8.
%
%   See also DPARAM, DPLOT, DINVMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: dmap.m 7 1998-05-10 04:37:19Z tad $

if isempty(zp)
  wp = [];
  return
end

n = length(z);
w = w(:);
beta = beta(:);
z = z(:);

% Quadrature data and error tolerance
if nargin < 6
  tol = 1e-8;
  qdat = scqdata(beta,8);
elseif length(qdat)==1
  tol = qdat;
  qdat = scqdata(beta,max(ceil(-log10(tol)),8));
else
  tol = 10^(-size(qdat,1));
end

shape = size(zp);
zp = zp(:);
zprow = zp.';
p = length(zp);
wp = zeros(p,1);

% For each point in zp, find nearest prevertex.
[dist,sing] = min(abs(zprow(ones(n,1),:) - z(:,ones(1,p))));
sing = sing(:);				% indices of prevertices

% Screen out images of prevertices
vertex = (dist(:) < tol);
wp(vertex) = w(sing(vertex));

% "Bad" points are closest to a prevertex of infinity.
atinf = find(isinf(w)); 		% infinite vertices
bad = ismember(sing,atinf) & ~vertex;

if any(bad)
  % Can't integrate starting at pre-infinity: find conformal center to use
  % as integration basis.
  if ~isinf(w(n-1))
    wc = w(n-1) + c*dquad(z(n-1),0,n-1,z,beta,qdat);
  else
    wc = w(n) + c*dquad(z(n),0,n,z,beta,qdat);
  end
end

% zs = the starting singularities
zs = z(sing);
% ws = map(zs)
ws = w(sing);

% Compute the map directly at "normal" points.
normal = ~bad & ~vertex;
if any(normal)
  I = dquad(zs(normal),zp(normal),sing(normal),z,beta,qdat);
  wp(normal) = ws(normal) + c*I;
end

% Compute map at "bad" points, using conformal center as basis, to avoid
% integration where right endpoint is too close to a singularity.
if any(bad)
  I = dquad(zp(bad),zeros(sum(bad),1),zeros(sum(bad),1),z,beta,qdat);
  wp(bad) = wc - c*I;
end

wp = reshape(wp,shape);
