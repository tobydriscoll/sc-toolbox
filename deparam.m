function [z,c,qdat] = deparam(w,beta,z0,options)
%DEPARAM Schwarz-Christoffel exterior parameter problem.  
%   [Z,C,QDAT] = DEPARAM(W,BETA) solves the Schwarz-Christoffel mapping
%   parameter problem with a disk as fundamental domain and the exterior
%   of the polygon specified by W as the target. W must be a vector of
%   the vertices of the polygon, specified in clockwise order, and BETA
%   should be a vector of the turning angles of the polygon; see
%   SCANGLES for details. If successful, DEPARAM will return Z, a vector
%   of the pre-images of W; C, the multiplicative constant of the
%   conformal map; and QDAT, an optional matrix of quadrature data used
%   by some of the other S-C routines.
%
%   [Z,C,QDAT] = DEPARAM(W,BETA,Z0) uses Z0 as an initial guess for Z.
%
%   [Z,C,QDAT] = DEPARAM(W,BETA,TOL) attempts to find an answer within
%   tolerance TOL. (Also see SCPAROPT.)
%       
%   [Z,C,QDAT] = DEPARAM(W,BETA,Z0,OPTIONS) uses a vector of control
%   parameters. See SCPAROPT.
%       
%   See also SCPAROPT, DRAWPOLY, DEDISP, DEPLOT, DEMAP, DEINVMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: deparam.m 129 2001-05-07 15:04:13Z driscoll $

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

err = sccheck('de',w,beta);
if err==1
  fprintf('Use SCFIX to make polygon obey requirements\n')
  error(' ')
end

[trace,tol,method] = parseopt(options);
if length(z0)==1
  tol = z0;
  z0 = [];
end
nqpts = max(ceil(-log10(tol)),2);
qdat = scqdata(beta,nqpts); 		% quadrature data
  
if n==2					% it's a slit
  z = [-1;1];

else
  % Set up normalized lengths for nonlinear equations
  len = abs(diff(w([n,1:n])));
  nmlen = abs(len(3:n-1)/len(2));
  
  % Set up initial guess
  if isempty(z0)
    y0 = zeros(n-1,1);
  else
    z0 = z0/z0(n);			% Fix z0(n)=1.
    th = angle(z0(:));
    th(th<=0) = th(th<=0) + 2*pi;
    dt = diff([0;th(1:n-1);2*pi]);
    y0 = log(dt(1:n-1)./dt(2:n));
  end
  
  % Solve nonlinear system of equations:

  % package data
  fdat = {n,beta,nmlen,qdat};
  % set options
  opt = zeros(16,1);
  opt(1) = trace;
  opt(2) = method;
  opt(6) = 100*(n-3);
  opt(8) = tol;
  opt(9) = min(eps^(2/3),tol/10);
  opt(12) = nqpts;
  try
    [y,termcode] = nesolve('depfun',y0,opt,fdat);
  catch
    % Have to delete the "waitbar" figure if interrupted
    close(findobj(allchild(0),'flat','Tag','TMWWaitbar'));
    error(lasterr)
  end
  if termcode~=1
    warning('Nonlinear equations solver did not terminate normally.')
  end
  
  % Convert y values to z
  cs = cumsum(cumprod([1;exp(-y)]));
  theta = 2*pi*cs/cs(n);
  z = ones(n,1);
  z(1:n-1) = [exp(i*theta(1:n-1))];
end

% Determine scaling constant
mid = z(1)*exp(i*angle(z(2)/z(1))/2);
c = (w(2) - w(1)) / (dequad(z(1),mid,1,z,beta,qdat)-...
    dequad(z(2),mid,2,z,beta,qdat));

