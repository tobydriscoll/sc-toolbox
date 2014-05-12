function [z,c,L,qdat] = rparam(w,beta,cnr,z0,options);
%RPARAM Schwarz-Christoffel rectangle parameter problem.
%   [Z,C,L,QDAT] = RPARAM(W,BETA,CORNERS) solves the Schwarz-Christoffel
%   parameter problem with a rectangle as fundamental domain and
%   interior of the specified polygon as the target. W must be a vector
%   of the vertices of the polygon, specified in counterclockwise
%   order. BETA is a vector of turning angles; see SCANGLES. CORNERS is
%   a 4-component vector specifying the indices of the vertices which
%   are the images of the corners of the rectangle. **Be sure** the
%   first two entries describe the LONG sides of the rectangle, and go
%   in counterclockwise order. If CORNERS is omitted, the user is
%   requested to select these vertices using the mouse.
%
%   If successful, RPARAM will return Z, a vector of the prevertices; C,
%   the multiplicative constant of the conformal map; L, a parameter
%   related to aspect ratio; and QDAT, an optional matrix of quadrature
%   data used by some of the other SC routines.
%
%   [Z,C,L,QDAT] = RPARAM(W,BETA,CORNERS,Z0) uses Z0 as an initial guess
%   for Z.  In this case, Z0 represents the image of prevertices on the
%   strip 0 <= Im z <= 1.  You can use R2STRIP to transform prevertices
%   from the rectangle to the strip.
%       
%   [Z,C,L,QDAT] = RPARAM(W,BETA,CORNERS,TOL) attempts to find an answer
%   within tolerance TOL. (Also see SCPAROPT.)
%
%   [Z,C,L,QDAT] = RPARAM(W,BETA,CORNERS,Z0,OPTIONS) uses a vector of
%   control parameters. See SCPAROPT.
%       
%   See also SCPAROPT, DRAWPOLY, RDISP, RPLOT, RMAP, RINVMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: rparam.m 298 2009-09-15 14:36:37Z driscoll $

n = length(w); 				% no. of vertices
w = w(:);
beta = beta(:);

% Set up defaults for missing args
if nargin < 5
  options = [];
  if nargin < 4
    z0 = [];
    if nargin < 3
      cnr = [];
    end
  end
end

[trace,tol,method] = parseopt(options);

if isempty(cnr)
  msg{1} = 'Select the images of the corners of the rectangle.';
  msg{2} = 'Go in counterclockwise order and select a long rectangle edge first.';
  cnr = scselect(w,beta,4,'Select corners',msg);
end

% Renumber vertices so that cnr(1)=1
renum = [cnr(1):n,1:cnr(1)-1];
w = w(renum);
beta = beta(renum);
cnr = rem(cnr-cnr(1)+1+n-1,n)+1;

if length(z0)==1
  tol = z0;
  z0 = [];
end
nqpts = max(ceil(-log10(tol)),4);
qdat = scqdata(beta,nqpts); 		% quadrature data

% Check input data.
err = sccheck('r',w,beta,cnr);
if err==1
  fprintf('Use SCFIX to make polygon obey requirements\n')
  error(' ')
end

atinf = (beta <= -1);

if isempty(z0)
  % Try to find a reasonable initial guess.
  dw = abs(diff(w([1:n,1])));		% side lengths
  dw(isinf(dw)) = mean(dw(~isinf(dw)))*ones(size(dw(isinf(dw))));
  % Estimate length and width (conformal modulus)
  len = mean([sum(dw(cnr(1):cnr(2)-1)), sum(dw(cnr(3):cnr(4)-1))]);
  wid = mean([sum(dw(cnr(2):cnr(3)-1)), sum(dw([cnr(4):n,1:cnr(1)-1]))]);
  modest = min(len/wid,100);
  % Evenly space prevertices to match this conformal modulus
  z0(cnr(1):cnr(2)) = linspace(0,modest,diff(cnr(1:2))+1);
  dx = z0(cnr(1)+1)-z0(cnr(1));
  z0(cnr(1)-1:-1:1) = z0(cnr(1))-dx*(1:cnr(1)-1);
  z0(cnr(2)+1:cnr(3)-1) = z0(cnr(2)) + dx*(1:diff(cnr(2:3))-1);
  z0(cnr(4):-1:cnr(3)) = i + linspace(0,modest,diff(cnr(3:4))+1);
  dx = z0(cnr(4)-1)-z0(cnr(4));
  z0(cnr(4)+1:n) = z0(cnr(4))-dx*(1:n-cnr(4));

else
  if length(z0)~=n
    error('Initial guess has wrong number of prevertices')
  end
  z0 = z0(renum);
  z0 = real(z0) + i*round(imag(z0));
  if any(imag(z0(cnr(1):cnr(2)))) | any(imag(z0(cnr(3):cnr(4)))==0)
    error('Initial guess has prevertices on wrong side of strip')
  end
end

% Convert z0 to unconstrained vars
y0 = zeros(n-3,1);
dz = diff(z0);
dz(cnr(3):n-1) = -dz(cnr(3):n-1);
y0(1:cnr(2)-2) = log(dz(1:cnr(2)-2));
y0(cnr(3)-1:cnr(4)-3) = log(dz(cnr(3)+1:cnr(4)-1));
y0(cnr(2)-1) = mean(log(dz([cnr(2)-1,cnr(3)])));

% Vertices on the "short" edges are transformed into the interval [-1,1],
% and then the Trefethen-style transformation is used.
L = z0(cnr(2)) - z0(cnr(1));
x = real(exp(pi*(L-conj(z0(cnr(2)+1:cnr(3)-1)))));
dx = -diff([1;x(:);-1]);
y0(cnr(2):cnr(3)-2) = log(dx(1:end-1)./dx(2:end));
x = real(exp(pi*z0(cnr(4)+1:n)));
dx = diff([-1;x(:);1]);
y0(cnr(4)-2:n-3) = log(dx(1:end-1)./dx(2:end));


% Find prevertices (solve param problem)

% Set up normalized lengths for nonlinear equations:
% indices of left and right integration endpoints
left = 1:n-2;				
right = 2:n-1;				
% delete indices corresponding to vertices at Inf
left(atinf(left)) = [];
right(atinf(right)) = [];
if atinf(n-1)
  right = [right,n];
end
cmplx = ((right-left) == 2);
% It's possible we replaced the last single condition by a complex one.
if length(cmplx)+sum(cmplx) > n-2
  cmplx(end) = 0;
end
% normalize lengths by w(2)-w(1)
nmlen = (w(right)-w(left))/(w(2)-w(1));
% abs value for finite ones
nmlen(~cmplx) = abs(nmlen(~cmplx));
% first entry is useless (=1)
nmlen(1) = [];

% Solve nonlinear system of equations:
% package data
%%figure
%%hp = line(NaN,NaN);
hp = [];
fdat = {n,beta,nmlen,left,right,logical(cmplx),qdat,cnr};
% set options
opt = zeros(16,1);
opt(1) = trace;
opt(2) = method;	
opt(6) = 100*(n-3);			% max # of iterations
opt(8) = tol;
opt(9) = min(eps^(2/3),tol/10);
opt(12) = nqpts;
try
  [y,termcode] = nesolve('rpfun',y0,opt,fdat);
catch
  % Have to delete the "waitbar" figure if interrupted
  close(findobj(allchild(0),'flat','Tag','TMWWaitbar'));
  error(lasterr)
end
if termcode~=1
  warning('Nonlinear equations solver did not terminate normally.')
end

% Convert y values to z on strip
z = rptrnsfm(y,cnr);
ends = find(diff(imag(z([1:n 1]))));
zs = [z(1:ends(1));Inf;z(ends(1)+1:ends(2));-Inf;z(ends(2)+1:n)];
bs = [beta(1:ends(1));0;beta(ends(1)+1:ends(2));0;beta(ends(2)+1:n)];
idx = [1:ends(1) n+1 ends(1)+1:ends(2) n+1 ends(2)+1:n n+1];
qs = qdat(:,[idx idx+n+1]);

% Determine multiplicative constant
mid = mean(zs(1:2));
g = stquad(zs(2),mid,2,zs,bs,qs) - stquad(zs(1),mid,1,zs,bs,qs);
c = -diff(w(1:2))/g;


% Find prevertices on the rectangle

% Find corners of rectangle
zs = z;
L = max( real(zs(cnr([2 3]))) ) - zs(cnr(1));
[K,Kp] = ellipkkp(L);
rect = [K;K+i*Kp;-K+i*Kp;-K];

l = logical(zeros(n,1));		% on left side 
l(cnr(3):cnr(4)) = 1;
r = logical(zeros(n,1));		% on right side
r(cnr(1):cnr(2)) = 1;
tl = (real(zs) > L) & imag(zs)>0;	% top-left side
tr = (real(zs) > L) & imag(zs)==0;	% top-right side
bl = (real(zs) < 0) & imag(zs)>0;	% bottom-left side
br = (real(zs) < 0) & imag(zs)==0;	% bottom-right side

% Initial guesses
z = zeros(size(zs));
% Corners
z(cnr) = rect;
% Left and right sides are simple
z(l) = linspace(rect(3),rect(4),diff(cnr(3:4))+1).';
z(r) = linspace(rect(1),rect(2),diff(cnr(1:2))+1).';
% Cluster the top and bottom guesses near 0
h = K/20;
z(tl) = i*Kp - h*(1:sum(tl));
z(tr) = i*Kp + h*(sum(tr):-1:1);
z(bl) = -h*(sum(bl):-1:1);
z(br) = h*(1:sum(br));

% Newton iteration on "r2strip() - zs = 0"
zn = z(:);
maxiter = 50;
done = logical(zeros(size(zn)));
done(cnr) = 1;
k = 0;  F = 0;
while ~all(done) & k < maxiter
  [F,dF] = r2strip(zn(~done),z(cnr),L);
  F = zs(~done) - F;

  % Adjust Newton step to stay exactly on rectangle boundary
  step = F./dF;
  lr = r(~done) | l(~done);			% on top or bottom
  step(lr) = i*imag(step(lr));
  step(~lr) = real(step(~lr));
  
  % Newton step
  znew = zn(~done) + step;
  
  % Keep prevertices from moving too far (past boundaries)
  % Left/right sides capped in Im direction
  x = min( max( imag(znew(lr)), 0 ), Kp);
  znew(lr) = real(znew(lr)) + i*x;
  % Top/bottom-left sides capped in Re direction
  tbl = tl(~done) | bl(~done);
  x = min( max(real(znew(tbl)), -K), -eps );
  znew(tbl) = i*imag(znew(tbl)) + x;
  % Top/bottom-right sides capped in Re direction
  tbr = tr(~done) | br(~done);
  x = min( max(real(znew(tbr)), eps), K );
  znew(tbr) = i*imag(znew(tbr)) + x;
  
  % Update
  zn(~done) = znew;
  done(~done) =  (abs(F) < tol);
  k = k + 1;
end
if any(abs(F)> tol)
  warning('Could not converge to the rectangle prevertices.')
end
z(:) = zn; 

% Undo renumbering
z(renum) = z;

