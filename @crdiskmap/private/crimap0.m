function zp = crimap0(wp,z,beta,aff,qdat,options)
%CRIMAP0 Single-embedding inverse map in crossratio formulation.
%   CRIMAP0(WP,Z,BETA,AFF) computes the inverse of WP under the map
%   defined by the prevertex embedding Z and the affine transformation
%   AFF(1:2).
%       
%   CRIMAP0(WP,Z,BETA,AFF,QDAT,OPTIONS) supplies the quadrature data
%   QDAT and the vector of OPTIONS as described in SCINVOPT.
%       
%   For more information on the algorithm, see DINVMAP. However, there
%   is no issue with starting points; the origin is always used as the
%   starting point.
%       
%   See also CRINVMAP, CRAFFINE, SCINVOPT, SCQDATA.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crimap0.m 298 2009-09-15 14:36:37Z driscoll $

% Parse input and initialize
import sctool.*
n = length(beta);
if nargin < 6
  options = [];
  if nargin < 5
    qdat = [];
  end
end
zp = 0*wp;
wp = wp(:);
lenwp = length(wp);

[ode,newton,tol,maxiter] = scinvopt(options);

if isempty(qdat)
  qdat = scqdata(beta,max(ceil(-log10(tol)),2));
end

% Ignore points with beta==0 in integration
mask = abs(beta) > eps;
beta2 = beta(mask);
z2 = z(mask);
n2 = sum(mask);

% Use origin of disk as initial guess for inverse images
w0 = aff(2)*ones(lenwp,1);
z0 = zeros(lenwp,1);
done = zeros(size(wp));

if ode
  % Use relaxed ODE tol if improving with Newton
  odetol = max(tol,1e-2*(newton));

  % Rescale dependent coordinate
  scale = (wp(~done) - w0(:));

  odefun = @(w,y) diskmap.imapfun(w,y,scale,z2,beta2,aff(1));
  [t,y] = ode23(odefun,[0,0.5,1],zeros(2*length(wp),1),odeset('abstol',odetol));
  [m,leny] = size(y);
  zp(:) = y(m,1:lenwp) + i*y(m,lenwp+1:leny);
  abszp = abs(zp);
  out = abszp > 1;
  zp(out) = zp(out)./abszp(out);
end
  
% Newton iterations
if newton
  % Setup
  if ~ode
    zn = z0(:);
    if length(z0)==1 & lenwp > 1
      zn = zn(:,ones(lenwp,1));
    end
    zn(done) = zp(done);
  else
    zn = zp(:);
  end

  % The iteration
  k = 0;
  while ~all(done) & k < 16
    F = wp(~done) - crmap0(zn(~done),z,beta,aff,qdat);
    m = length(F);
    dF = aff(1)*exp(sum(beta2(:,ones(m,1)).*...
	log(1-(zn(~done,ones(n2,1)).')./z2(:,ones(m,1)))));
    zn(~done) = zn(~done) + F(:)./dF(:);
    done(~done) = (abs(F)< tol);
    k = k+1;
  end
  if any(abs(F) > tol)
    disp('Warning in crinvmap: Solution may be inaccurate')
    fprintf('Maximum residual = %.3g\n',max(abs(F)))
  end
  zp(:) = zn; 
end

end
