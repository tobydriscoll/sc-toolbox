function [z,zb,c,qdat] = rsparam(w,beta,branch,z0,options)
%RSPARAM Schwarz-Christoffel Riemann surface parameter problem.

%   [Z,C,QDAT] = RSPARAM(W,BETA) solves the Schwarz-Christoffel mapping
%   parameter problem with the disk as fundamental domain and the
%   polygon specified by W as the target. W must be a vector of the
%   vertices of the polygon, specified in counterclockwise order, and
%   BETA should be a vector of the turning angles of the polygon; see
%   SCANGLE for details. The polygon may be multiply sheeted, in which
%   case the branch points must be given.

%   If successful, DPARAM will return Z, a vector
%   of the pre-images of W; C, the multiplicative constant of the
%   conformal map; and QDAT, an optional matrix of quadrature data used
%   by some of the other S-C routines.
%
%   [Z,C,QDAT] = RSPARAM(W,BETA,Z0) uses Z0 as an initial guess for Z.
%       
%   [Z,C,QDAT] = RSPARAM(W,BETA,TOL) attempts to find an answer within
%   tolerance TOL. (Also see SCPAROPT.)
%
%   [Z,C,QDAT] = RSPARAM(W,BETA,Z0,OPTIONS) uses a vector of control
%   parameters. See SCPAROPT.



%   Copyright 2002 by Toby Driscoll.
%   $Id: rsparam.m 298 2009-09-15 14:36:37Z driscoll $

n = length(w);				% no. of vertices
w = w(:);
beta = beta(:);
B = length(branch);

% Set up defaults for missing args
if nargin < 5
  options = [];
  if nargin < 4
    z0 = [];
  end
end

err = sccheck('d',w,beta);
if err==1
  fprintf('Use SCFIX to make polygon obey requirements\n')
  error(' ')
end

[trace,tol,method] = scparopt(options);
if length(z0)==1
  tol = z0;
  z0 = [];
end
nqpts = max(ceil(-log10(tol)),4);
qdat = scqdata(beta,nqpts); 		% quadrature data

atinf = (beta <= -1);

if n==3
  % Trivial solution
  z = [-1i;(1-1i)/sqrt(2);1];

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
%    y0 = randn(n-3+2*B,1);%zeros(n-3+2*B,1);%
    y0=[zeros(n-3,1);randn(2*B,1)];
  else
    zb = z0{2};
    z0 = z0{1};
    z0 = z0(:)./abs(z0(:));
    % Moebius to make th(n-2:n)=[1,1.5,2]*pi;
    Am = moebius(z0(n-2:n),[-1;-1i;1]);
    z0 = Am(z0);
    th = angle(z0);
    th(th<=0) = th(th<=0) + 2*pi;
    dt = diff([0;th(1:n-2)]);
    y0 = log(dt(1:n-3)./dt(2:n-2));
    % Map branch images to upper HP
    Am = moebius(1i,1i,1,-1);
    u = Am(zb);
    u0 = [real(u) sqrt(imag(u))]';
%%    u = sqrt(Am(zb));
%%    % Transform positive numbers using log.
%%    u0 = [log(real(u));log(imag(u))];
    y0 = [y0;u0(:)];
  end
  
  % Solve nonlinear system of equations:
  
  % package data
  fdat = {w,beta,nmlen,branch,left,right,logical(cmplx),qdat};
  % set options
  opt = zeros(16,1);
  opt(1) = trace;
  opt(2) = method;
  opt(6) = 100*(n-3);
  opt(8) = tol;
  opt(9) = min(eps^(2/3),tol/10);
  opt(12) = nqpts;
  [~,termcode] = nesolve(@rspfun,y0,opt,fdat);
  if termcode~=1
    disp('Warning: Nonlinear equations solver did not terminate normally')
  end
  y = fsolve(@rspfun,y0,...
             optimset('display','iter','tolfun',tol,'large','off','tolx',tol),fdat);

  % Convert y values to z
  [z,zb] = vartsfm(y,n);
end

% Determine scaling constant
mid = (z(1)+z(2))/2;
c = (w(1) - w(2))/...
    (rsquad(z(2),mid,2,z,beta,zb,qdat) - rsquad(z(1),mid,1,z,beta,zb,qdat));


function F = rspfun(y,fdat)
%   Returns residual for solution of nonlinear equations.

[w,beta,nmlen,branch,left,right,cmplx,qdat] = deal(fdat{:});
n = length(w);

% Convert y values to z (prevertices)
[z,zb] = vartsfm(y,n);
theta = angle(z(1:n-3));

% Check crowding.
if any(diff(theta)<eps) || any(isnan(theta))
  % Since abs(y) is large, use it as the penalty function.
  F = y;
  disp('Warning: Severe crowding')
  return
end

% Compute the integrals appearing in nonlinear eqns.
zleft = z(left);
zright = z(right);
angl = angle(zleft);
mid = exp(1i*(angl + rem(angle(zright./zleft)+2*pi,2*pi)/2));
% For integrals between nonadjacent singularities, choose 0 as intermediate
% integration point.
mid(cmplx) = zeros(size(mid(cmplx)));
% If any are complex, the first one must be too.
if any(cmplx)
  cmplx(1) = 1;
end

ints = rsquad(zleft,mid,left,z,beta,zb,qdat) - ...
       rsquad(zright,mid,right,z,beta,zb,qdat);

%%if any(cmplx)
%%  ints(cmplx) = rsquad(zleft(cmplx),mid(cmplx),left(cmplx),...
%%                       z,beta,zb,qdat) - ...
%%      rsquad(zright(cmplx),mid(cmplx),right(cmplx),z,beta,zb,qdat);
%%end

% Branch point images
dist = abs(z(:,ones(length(zb),1))-zb(:,ones(n,1)).');
dist(isinf(w),:) = inf;
[m,idx] = min(dist,[],1);
Q = rsquad(z(idx),zb,idx,z,beta,zb,qdat);


if any(ints==0)
  % Singularities were too crowded in practice.
  F = y;
  disp('Warning: Severe crowding')
else
  % Compute nonlinear equation residual values.
  cmplx(1) = 0;
  F1 = ints(~cmplx); 		
  F1 = abs(F1(2:end)/F1(1));
  F2 = ints(cmplx)/ints(1);
  F = [F1;real(F2);imag(F2)] - nmlen;
  G = (branch-w(idx))/(w(right(1))-w(left(1))) - Q/ints(1);
  F = [F;real(G);imag(G)];
%%  G = [(branch-w(idx))/(w(right(1))-w(left(1)))  Q/ints(1)];
%%  F = [F;log([real(G(:,1))./real(G(:,2));imag(G(:,1))./imag(G(:,2))])];
end


function [z,zb] = vartsfm(y,n)

cs = cumsum(cumprod([1;exp(-y(1:n-3))]));
theta = pi*cs(1:n-3)/cs(length(cs));
z = ones(n,1);
z(1:n-3) = exp(1i*theta);
z(n-2:n-1) = [-1;-1i];

u = y(n-2:2:end);
v = y(n-1:2:end).^2;
zb = (u+1i*v-1i)./(u+1i*v+1i);
