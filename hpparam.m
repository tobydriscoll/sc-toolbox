function [z,c,qdat] = hpparam(w,beta,z0,options);
%HPPARAM Schwarz-Christoffel half-plane parameter problem.
%   [Z,C,QDAT] = HPPARAM(W,BETA) solves the Schwarz-Christoffel
%   parameter problem with the upper half-plane as fundamental domain
%   and interior of the specified polygon as the target. W must be a
%   vector of the vertices of the polygon, specified in counterclockwise
%   order. BETA is a vector of turning angles; see SCANGLES. If
%   successful, HPPARAM will return Z, a vector of the pre-images of W;
%   C, the multiplicative constant of the conformal map; and QDAT, an
%   optional matrix of quadrature data used by some of the other S-C
%   routines.
%
%   [Z,C,QDAT] = HPPARAM(W,BETA,Z0) uses Z0 as an initial guess for Z.
%       
%   [Z,C,QDAT] = HPPARAM(W,BETA,TOL) attempts to find an answer within
%   tolerance TOL. (Also see SCPAROPT.)
%
%   [Z,C,QDAT] = HPPARAM(W,BETA,Z0,OPTIONS) uses a vector of control
%   parameters. See SCPAROPT.
%       
%   See also SCPAROPT, DRAWPOLY, HPDISP, HPPLOT, HPMAP, HPINVMAP.

%   Copyright 1998--2001 by Toby Driscoll.
%   $Id: hpparam.m 298 2009-09-15 14:36:37Z driscoll $

n = length(w); 				% no. of vertices
w = w(:);
beta = beta(:);

% Set up defaults for missing args
if nargin < 4
  options = [];
  if nargin < 3
    z0 = [];
  end
end

err = sccheck('hp',w,beta);
if err==1
  fprintf('Use SCFIX to make polygon obey requirements\n')
  error(' ')
end

[trace,tol,method] = parseopt(options);
if length(z0)==1
  tol = z0;
  z0 = [];
end
nqpts = max(ceil(-log10(tol)),4);
qdat = scqdata(beta(1:n-1),nqpts); 	% quadrature data

atinf = (beta <= -1);

% Find prevertices (solve param problem)
if n==3
  z = [-1;1;Inf];

else

  % Set up normalized lengths for nonlinear equations:

  % indices of left and right integration endpoints
  left = find(~atinf(1:n-2));				
  right = 1+find(~atinf(2:n-1));				
  cmplx = ((right-left) == 2);
  % normalize lengths by w(2)-w(1)
  nmlen = (w(right)-w(left))/(w(2)-w(1));
  % abs value for finite ones; Re/Im for infinite ones
  nmlen = [abs(nmlen(~cmplx));real(nmlen(cmplx));imag(nmlen(cmplx))];
  % first entry is useless (=1)
  nmlen(1) = [];
  
  % Set up initial guess
  if isempty(z0)
    z0 = linspace(-1,1,n-1)';
  else
    z0 = z0(:);
    z0 = z0*2/(z0(n-1)-z0(1));
    z0 = z0-z0(1)-1;
  end
  %y0 = log(diff(z0(1:n-2)));
  y0 = log(diff(z0(1:n-2))./diff(z0(2:n-1)));

  % Solve nonlinear system of equations:

  % package data
  fdat = {n,beta(1:n-1),nmlen,left,right,logical(cmplx),qdat};
  % set options
  opt = zeros(16,1);
  opt(1) = trace;
  opt(2) = method;
  opt(6) = 100*(n-3);
  opt(8) = tol;
  opt(9) = min(eps^(2/3),tol/10);
  opt(12) = nqpts;
  try
    [y,termcode] = nesolve('hppfun',y0,opt,fdat);
  catch
    % Have to delete the "waitbar" figure if interrupted
    close(findobj(allchild(0),'flat','Tag','TMWWaitbar'));
    error(lasterr)
  end
  if termcode~=1
    warning('Nonlinear equations solver did not terminate normally.')
  end

  % Convert y values to z
  %z = [cumsum([-1;exp(y)]);1;Inf];
  cp = cumprod([1;exp(-y)]);
  z = [0;cumsum(cp)] - [flipud(cumsum(flipud(cp)));0];
  z = [z/z(n-1);Inf];
end

% Determine multiplicative constant
mid = mean(z(1:2));
g = hpquad(z(2),mid,2,z(1:n-1),beta(1:n-1),qdat) -...
    hpquad(z(1),mid,1,z(1:n-1),beta(1:n-1),qdat);
c = (w(1)-w(2))/g;


function F = hppfun(y,fdat)
%   Returns residual for solution of nonlinear equations. 

n = fdat(1,1);
beta = fdat(1:n-1,2);
nmlen = fdat(1:n-3,3);
rows = 1:fdat(2,1);
left = fdat(rows,4);
right = fdat(rows,5);
cmplx = logical(fdat(rows,6));
qdat = fdat(1:fdat(3,1),7:fdat(4,1));

% Transform y (unconstr. vars) to z (prevertices)
%z = [cumsum([-1;exp(y)]);1];
cp = cumprod([1;exp(-y)]);
z = [0;cumsum(cp)] - [flipud(cumsum(flipud(cp)));0];
z = z/z(n-1);

% Check crowding of singularities.
if any(diff(z)<eps) | any(isinf(z))
  % Since abs(y) is large, use it as the penalty function.
  F = y;
  disp('Warning: Severe crowding')
  return
end

% Compute the integrals appearing in nonlinear eqns.
zleft = z(left);
zright = z(right);
mid = mean([zleft.' ; zright.']).';
% For integrals between non-adjacent singularities, choose intermediate
% points in the upper half-plane.
mid(cmplx) = mid(cmplx) + i*(zright(cmplx)-zleft(cmplx))/2;
ints = hpquad(zleft,mid,left,z,beta,qdat) - ...
    hpquad(zright,mid,right,z,beta,qdat);

if any(ints==0)
  % Singularities were too crowded in practice.
  F = y;
  disp('Warning: Severe crowding')
else
  % Compute nonlinear equation residual values.
  F1 = abs(ints(~cmplx));		% F1(1) = abs(ints(1))
  F1 = F1(2:length(F1))/F1(1);
  F2 = ints(cmplx)/ints(1);
  F = [F1;real(F2);imag(F2)] - nmlen;
end
