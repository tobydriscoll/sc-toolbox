function wp = hpmap(zp,w,beta,z,c,qdat)
%HPMAP  Schwarz-Christoffel half-plane map.
%   HPMAP(ZP,W,BETA,Z,C,QDAT) computes the values of the
%   Schwarz-Christoffel half-plane map at the points in vector ZP.  The
%   polygon's vertices should be given in W and the arguments Z, C, and
%   QDAT should be computed by HPPARAM.  HPMAP returns a vector the same
%   size as ZP.
%   
%   HPMAP(ZP,W,BETA,Z,C,TOL) uses quadrature data intended to give an
%   answer accurate to within TOL.
%   
%   HPMAP(ZP,W,BETA,Z,C) uses a tolerance of 1e-8.
%
%   See also HPPARAM, HPPLOT, HPINVMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: hpmap.m 298 2009-09-15 14:36:37Z driscoll $

if isempty(zp)
  wp = [];
  return
end
import sctool.*

n = length(w);
w = w(:);
beta = beta(:);
z = z(:);

% Quadrature data and error tolerance
if nargin < 6
  tol = 1e-8;
  qdat = scqdata(beta(1:n-1),8);
elseif length(qdat)==1
  tol = qdat;
  qdat = scqdata(beta(1:n-1),max(ceil(-log10(tol)),8));
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
zpinf = isinf(zp);
wp(zpinf) = w(n);
vertex = vertex | zpinf;

% "Bad" points are closest to a prevertex of infinity.
atinf = find(isinf(w)); 		% infinite vertices
bad = ismember(sing,atinf) & ~vertex;

if any(bad)
  % Can't integrate starting at pre-infinity: which neighboring prevertex
  % to use?
  direcn = real(zp(bad)-z(sing(bad)));
  sing(bad) = sing(bad) + sign(direcn) + (direcn==0);
  % Midpoints of these integrations 
  mid = (z(sing(bad)) + zp(bad)) / 2;
end
  
% zs = the starting singularities
zs = z(sing);
% ws = f(zs)
ws = w(sing);

% Compute the map directly at "normal" points.
normal = ~bad & ~vertex;
if any(normal)
  I = hpquad(zs(normal),zp(normal),sing(normal),z(1:n-1),beta(1:n-1),qdat);
  wp(normal) = ws(normal) + c*I;
end

% Compute map at "bad" points, in stages. Stop at midpoint to avoid
% integration where right endpoint is close to a singularity.
if any(bad)
  I1 = hpquad(zs(bad),mid,sing(bad),z(1:n-1),beta(1:n-1),qdat);
  I2 = -hpquad(zp(bad),mid,zeros(sum(bad),1),z(1:n-1),beta(1:n-1),qdat);
  wp(bad) = ws(bad) + c*(I1 + I2);
end

wp = reshape(wp,shape);
