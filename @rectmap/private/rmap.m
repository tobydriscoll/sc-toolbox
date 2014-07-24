function wp = rmap(zp,w,beta,z,c,L,qdat)
%RMAP   Schwarz-Christoffel rectangle map.
%   RMAP(ZP,W,BETA,Z,C,L,QDAT) computes the values of the
%   Schwarz-Christoffel rectangle map at the points in vector ZP.  The
%   remaining arguments are as in RPARAM.  RMAP returns a vector the
%   same size as ZP.
%   
%   RMAP(ZP,W,BETA,Z,C,L,TOL) uses quadrature data intended to give an
%   answer accurate to within TOL.
%   
%   RMAP(ZP,W,BETA,Z,C,L) uses a tolerance of 1e-8.
%
%   See also RPARAM, RPLOT, RINVMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: rmap.m 298 2009-09-15 14:36:37Z driscoll $

if isempty(zp)
  wp = [];
  return
end

import sctool.*
n = length(w);
wp = z;
w = w(:);
beta = beta(:);
z = z(:);
[w,beta,z,corners] = rcorners(w,beta,z);

if nargin < 7
  qdat = scqdata(beta,8);
elseif length(qdat)==1
  qdat = scqdata(beta,max(ceil(-log10(qdat)),8));
end

% Map prevertices to strip
K = max(real(z));
Kp = max(imag(z));
zs = r2strip(z,z,L);
zs = real(zs) + i*round(imag(zs));	% put them *exactly* on edges

% Add in ends of strip
ends = find(diff(imag(zs([1:n 1]))));
zs = [zs(1:ends(1));Inf;zs(ends(1)+1:ends(2));-Inf;zs(ends(2)+1:n)];
ws = [w(1:ends(1));NaN;w(ends(1)+1:ends(2));NaN;w(ends(2)+1:n)];
bs = [beta(1:ends(1));0;beta(ends(1)+1:ends(2));0;beta(ends(2)+1:n)];
% Extend qdat with useless columns at ends
idx = [1:ends(1) n+1 ends(1)+1:ends(2) n+1 ends(2)+1:n n+1];
qdat = qdat(:,[idx idx+n+1]);

wp = zeros(size(zp));
zp = zp(:);
p = length(zp);

% Trap points which map to +/-Inf on the strip.
bad = abs(zp) < 2*eps;
zp(bad) = zp(bad) + 100*eps;
bad = abs(zp-i*Kp) < 2*eps;
zp(bad) = zp(bad) - i*100*eps*Kp;

% Map from rectangle to strip.
yp = r2strip(zp,z,L);

% Now map from strip to polygon.
wp(:) = stripmap.evaluate(yp,ws,bs,zs,c,qdat);


