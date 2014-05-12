function [z,c,qdat] = stparam(w,beta,ends,z0,options);
%STPARAM Schwarz-Christoffel strip parameter problem.
%   [Z,C,QDAT] = STPARAM(W,BETA,ENDS) solves the Schwarz-Christoffel
%   parameter problem with the infinite strip as fundamental domain and
%   interior of the specified polygon as the target. W must be a vector
%   of the vertices of the polygon, specified in counterclockwise
%   order. BETA is a vector of turning angles; see SCANGLES. ENDS is a
%   2-vector whose entries are the indices of the vertices which are the
%   images of the left and right ends of the strip. If ENDS is omitted,
%   the user is requested to select these vertices using the mouse.
%
%   If successful, STPARAM will return Z, a vector of the pre-images of
%   W; C, the multiplicative constant of the conformal map; and QDAT, an
%   optional matrix of quadrature data required by some of the other SC
%   routines.
%
%   [Z,C,QDAT] = STPARAM(W,BETA,ENDS,Z0) uses Z0 as an initial guess for
%   Z.
%       
%   [Z,C,QDAT] = STPARAM(W,BETA,ENDS,TOL) attempts to find an answer
%   within tolerance TOL. (Also see SCPAROPT.)
%
%   [Z,C,QDAT] = STPARAM(W,BETA,ENDS,Z0,OPTIONS) uses a vector of
%   control parameters.  See SCPAROPT.
%       
%   See also SCPAROPT, DRAWPOLY, STDISP, STPLOT, STMAP, STINVMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: stparam.m 298 2009-09-15 14:36:37Z driscoll $

% Set up defaults for missing args
if nargin < 5
  options = [];
  if nargin < 4
    z0 = [];
    if nargin < 3
      ends = [];
    end
  end
end

[trace,tol,method] = parseopt(options);

if isempty(ends)
  msg = 'Select the vertices that map to the ends of the strip.';
  ends = scselect(w,beta,2,'Select ends',msg);
end

N = length(w); 				% no. of vertices
w = w(:);
beta = beta(:);
% Renumber vertices so that the ends of the strip map to w([1,k])
renum = [ends(1):N,1:ends(1)-1];
w = w(renum);
beta = beta(renum);
k = find(renum==ends(2));
% n: number of finite prevertices
n = N-2;
% nb: number of prevertices on bottom edge of strip
nb = k-2;

% Check input data.
err = sccheck('st',w,beta,ends);
if err==1
  fprintf('Use SCFIX to make polygon obey requirements\n')
  error(' ')
end

if length(z0)==1
  tol = z0;
  z0 = [];
end
nqpts = max(ceil(-log10(tol)),4);
qdat = scqdata(beta,nqpts); 		% quadrature data
atinf = (beta <= -1);

% Ignore images of ends of strip.
w([1,k]) = [];
atinf([1,k]) = [];

if isempty(z0)
  % Make initial guess based on polygon.
  z0 = zeros(n,1);
  if any(atinf)
    % Can't base it on relative side lengths.
    scale = (abs(w(nb)-w(1))+abs(w(n)-w(nb+1)))/2;
    z0(1:nb) = linspace(0,scale,nb)';
    z0(nb+1:n) = i + flipud(linspace(0,scale,n-nb)');
  else
    % This is from Louis Howell's code.
    scale = (abs(w(n)-w(1))+abs(w(nb)-w(nb+1)))/2;
    z0(1:nb) = cumsum([0;abs(w(2:nb)-w(1:nb-1))]/scale);
    if nb+1==n
      z0(n) = mean(z0([1,nb]));
    else
      z0(n:-1:nb+1) = cumsum([0;abs(w(n:-1:nb+2)-w(n-1:-1:nb+1))]/scale);
    end
    scale = sqrt(z0(nb)/z0(nb+1));
    z0(1:nb) = z0(1:nb)/scale;
    z0(nb+1:n) = i + z0(nb+1:n)*scale;
  end
else
  z0 = z0(renum);
  if length(z0)==N 
    if ~all(isinf(z0([1,k])))
      error('Starting guess does not match ends of strip')
    end
    z0([1,k]) = [];
  elseif length(z0)==n-1
    z0 = [0;z0];
  end
end
y0 = [log(diff(z0(1:nb)));real(z0(nb+1));log(-diff(z0(nb+1:n)))];

% Find prevertices (solve param problem)

% Set up normalized lengths for nonlinear equations:
% indices of left and right integration endpoints
left = [1,1:n-1]';				
right = [n,2:n]';				
% delete indices corresponding to vertices at Inf
left([find(atinf)+1;nb+1]) = [];
right([find(atinf);nb+1]) = [];
cmplx = ((right-left) == 2);
cmplx(1) = 0;
cmplx(2) = 1;
% normalize lengths 
nmlen = (w(right)-w(left))/(w(n)-w(1));
% abs value for finite ones
nmlen(~cmplx) = abs(nmlen(~cmplx));
% first entry is useless (=1)
nmlen(1) = [];
  
% Solve nonlinear system of equations:

% package data
fdat = {n,nb,beta,nmlen,left,right,logical(cmplx),qdat};
% set options
opt = zeros(16,1);
opt(1) = trace;
opt(2) = method;
opt(6) = 100*(n-1);
opt(8) = tol;
opt(9) = tol/10;
opt(12) = nqpts;
try
  [y,termcode] = nesolve('stpfun',y0,opt,fdat);
catch
  % Have to delete the "waitbar" figure if interrupted
  close(findobj(allchild(0),'flat','Tag','TMWWaitbar'));
  error(lasterr)
end
if termcode~=1
  warning('Nonlinear equations solver did not terminate normally.')
end


% Convert y values to z
z = zeros(n,1);
z(2:nb) = cumsum(exp(y(1:nb-1)));
z(nb+1:n) = i+cumsum([y(nb);-exp(y(nb+1:n-1))]);
z = [-Inf;z(1:nb);Inf;z(nb+1:n)];

% Determine multiplicative constant
mid = mean(z(2:3));
g = stquad(z(3),mid,3,z,beta,qdat) - stquad(z(2),mid,2,z,beta,qdat);
c = (w(1)-w(2))/g;

% Undo renumbering
z(renum) = z;
qdat(:,renum) = qdat(:,1:N);
qdat(:,N+1+renum) = qdat(:,N+1+(1:N));
