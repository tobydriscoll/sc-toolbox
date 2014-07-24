function [zp,flag] = dinvmap(wp,w,beta,z,c,qdat,z0,options)
%DINVMAP Schwarz-Christoffel disk inverse map.
%   DINVMAP(WP,W,BETA,Z,C,TOL) computes the inverse of the
%   Schwarz-Christoffel disk map (i.e., from a polygon to the disk) at
%   the points given in vector WP. The other arguments are as in
%   DPARAM. TOL is a scalar tolerance, or a quadrature-data matrix QDAT
%   as returned by SCQDATA, or may be omitted.
%       
%   The default algorithm is to solve an ODE in order to obtain a fair
%   approximation for ZP, and then improve ZP with Newton iterations.
%   The ODE solution at WP requires a vector Z0 whose forward image W0
%   is such that for each j, the line segment connecting WP(j) and W0(j)
%   lies inside the polygon.  By default Z0 is chosen by a fairly robust
%   automatic process.  Using a parameter (see below), you can choose to
%   use either an ODE solution or Newton iterations exclusively.
%
%   DINVMAP(WP,W,BETA,Z,C,TOL,Z0) has two interpretations.  If the ODE
%   solution is being used, Z0 overrides the automatic selection of
%   initial points.  (This can be handy in convex polygons, where the
%   choice of Z0 is trivial.)  Otherwise, Z0 is taken as an initial
%   guess to ZP.  In either case, if length(Z0)==1, the value Z0 is used
%   for all elements of WP; otherwise, length(Z0) should equal
%   length(WP).
%
%   DINVMAP(WP,W,BETA,Z,C,TOL,Z0,OPTIONS) uses a vector of parameters
%   that control the algorithm.  See SCINVOPT.
%
%   [ZP,FLAG] = DINVMAP(...) also returns a vector of indices where the
%   method was unable to produce a sufficiently small residual. A warning
%   is issued when this occurs.
% 
%   See also SCINVOPT, DPARAM, DMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: dinvmap.m 279 2007-05-14 20:14:24Z driscoll $

n = length(w);
beta = beta(:);
z = z(:);
zp = zeros(size(wp));
wp = wp(:);
lenwp = length(wp);

if nargin < 8
  options = [];
  if nargin < 7
    z0 = [];
    if nargin < 6
      qdat = [];
    end
  end
end

[ode,newton,tol,maxiter] = sctool.scinvopt(options);

if isempty(qdat)
  qdat = tol;
end

if length(qdat)==1
  qdat = sctool.scqdata(beta,max(ceil(-log10(qdat)),2));
end

done = zeros(size(wp));
% First, trap all points indistinguishable from vertices, or they will cause
% trouble.
% Modified 05/14/2007 to work around bug in matlab 2007a.
for j=1:n
  idx = find(abs(wp-w(j)) < 3*eps);
  zp(idx) = z(j);
  done(idx) = 1;
end
lenwp = lenwp - sum(done);
if lenwp==0, flag = []; return, end

% ODE
if ode
  if isempty(z0)
    % Pick a value z0 (not a singularity) and compute the map there.
    map = @(zp) dmap(zp,w,beta,z,c,qdat);
    [z0,w0] = sctool.findz0('d',wp(~done),map,w,beta,z,c,qdat);
  else
    w0 = dmap(z0,w,beta,z,c,qdat);
    if length(z0)==1 && lenwp > 1
      z0 = z0(:,ones(lenwp,1)).';
      w0 = w0(:,ones(lenwp,1)).';
    end
    w0 = w0(~done);
    z0 = z0(~done);
  end

  % Use relaxed ODE tol if improving with Newton.
  odetol = max(tol,1e-4*(newton));
  opt = odeset('abstol',odetol,'reltol',odetol);

  % Rescale dependent coordinate
  scale = (wp(~done) - w0(:));

  % Solve ODE
  z0 = [real(z0);imag(z0)];
  odefun = @(w,y) dimapfun(w,y,scale,z,beta,c);
  [t,y] = ode113(odefun,[0,0.5,1],z0,opt);
  [m,leny] = size(y);
  zp(~done) = y(m,1:lenwp)+sqrt(-1)*y(m,lenwp+1:leny);
  abszp = abs(zp);
  out = abszp > 1;
  zp(out) = zp(out)./abszp(out);
end

% Newton iterations
if newton
  if ~ode
    zn = z0(:);
    if length(z0)==1 && lenwp > 1
      zn = zn(:,ones(lenwp,1));
    end
    zn(done) = zp(done);
  else
    zn = zp(:);
  end
    
  wp = wp(:);
  k = 0;
  while ~all(done) && k < maxiter
    F = wp(~done) - dmap(zn(~done),w,beta,z,c,qdat);
    m = length(F);
    dF = c*exp(sum(beta(:,ones(m,1)).*...
          log(1-(zn(~done,ones(n,1)).')./z(:,ones(m,1)))));
    zn(~done) = zn(~done) + F(:)./dF(:);
    out = abs(zn) > 1;
    zn(out) = sign(zn(out));
    done(~done) = (abs(F)< tol);
    k = k+1;
  end
  if any(abs(F)> tol)
    str = sprintf('Check solution; maximum residual = %.3g\n',max(abs(F)));
    warning(str)
  end
  zp(:) = zn; 
end;

flag = find(~done);

end

