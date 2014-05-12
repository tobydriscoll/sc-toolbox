function [phi,fd,fr,sidenum] = lapsolve(p,bdata)
%LAPSOLVE Solve Laplace's equation on a polygon.
%
%   PHI = LAPSOLVE(P,BDATA) produces a solution to Laplace's equation
%   with piecewise constant boundary conditions on a polygon P. The
%   vector BDATA should have one entry for each side of P: A constant
%   value indicates a Dirichlet condition, while NaN indicates a
%   homogeneous Neumann condiiton.
%
%   PHI = LAPSOLVE(P) opens a GUI for selecting the boundary
%   conditions interactively.
%
%   P may instead be an HPLMAP, in which case the target polygon of
%   the map is used. This can save effort if multiple BVPs are to be
%   solved on the same region.
%
%   The output argument is a COMPOSITE function. The syntax PHI(z)
%   should be used to compute PHI at points. It cannot be inverted or
%   differentiated explicitly. To plot the solution you can use
%
%      phi = lapsolve(p);
%      [tri,x,y] = triangulate(p);
%      trisurf(tri,x,y,phi(x+i*y))
%
%   [PHI,FD,FR,SIDENUM] = LAPSOLVE(P,...) returns additional
%   information. FD and FR are two Schwarz-Christoffel maps (of type
%   DISKMAP and RIESURFMAP) such that
%
%                           -1
%      phi(z) = real( fr( fd  (z) ) ).
%
%   The polygon Q=POLYGON(FR) is conformally equivalent to P but may
%   have additional trivial vertices. SIDENUM helps with the
%   identification. SIDENUM(k) is the shows which side of P that side
%   k of Q lies on. For example, BDATA(SIDENUM) extends the boundary
%   data to Q.
%
%   See also DISKMAP, RIESURFMAP, COMPOSITE, POLYGON/TRIANGULATE.


% First, map to half-plane if not done.
if isequal(class(p),'hplmap')
  f = p;
  p = polygon(f);
else
  err = sccheck('hp',vertex(p),angle(p)-1);
  if err
    fprintf(['Polygon does not meet HPLMAP requirements. Use SCFIX first' ...
	     ' and adjust boundary data accordingly.']);
    error(' ');
  end
  f = hplmap(p);
end

% Get BC if not specified.
if nargin==1
  bdata = lapsolvegui(p);
end

% Interpret Neumann and Dirichlet selections.
bdata = bdata(:);
n = length(bdata);
if n~=length(p)
  error('Incorrect number of boundary conditions.')
end
neumann = isnan(bdata);
dirichlet = ~neumann;

% === Linear parameter problem ===

% Prevertices. 
z = prevertex(f);
% Translate them to avoid zero.
z = z-3;

% Angles of target.
d = diff(dirichlet([n 1:n]));
alpha = 0.5 * (d~=0);     % D/N or N/D transition
NN = (d==0) & neumann;    % adjacent N sides
alpha(NN) = 1;            % no turn
DD = (d==0) & (bdata([n 1:n-1])==bdata);  % adj D sides with same value
alpha(DD) = 1;
dirichlet = dirichlet & ~DD;
diridx = find(dirichlet);

beta = alpha - 1;

% Index of problem.
kappa = -sum( beta );
if kappa < 2
  error('At least two distinct Dirichlet values must be prescribed.')
end

% Endpoints of integration (points on the Dirichlet sides).
Z = [ z(1:n-1); 0; -10];
ze = 0.5*( Z(diridx) + Z(diridx+1) );

% Global rotation.
if neumann(n-1) 
  rho = 1;
else 
  rho = i;
end

% Integrate between neighboring midpoints. 
zs = [ z(1:n-1); 0];  % includes generic monomial
bs = [ beta(1:n-1); 0];
M = zeros(kappa-1,kappa-1);
s = zeros(size(ze));  % no singular endpoints.
for d = 1:kappa-1
  qdata = scqdata(bs,15);
  M(:,d) = hpquad(ze(1:end-1),ze(2:end),s,s,zs,bs,qdata);
  bs(end) = bs(end) + 1;                  % next monomial
end

% Find the slit polynomial (high degree last).
dphi = diff( bdata(dirichlet) );
poly = real(rho*M) \ dphi;

% Slit and branch locations.
r = roots(flipud(poly));
% Roots in lower half-plane are redundant.
r(imag(r)<-1e-12) = [];
idx = imag(r) > 1e-12;
zb = r(idx);
slit = real(r(~idx));

% Final prevertices and angles.
ns = length(slit);
z = [ z; slit; ];
alpha = [ alpha; 2*ones(ns,1) ];
% Keep track of which original side was the "parent", so we know what BC
% is applied where. (To be filled in after sorting.)
sidenum = [ (1:n)'; NaN*ones(ns,1) ];       

% Sort.
[z,idx] = sort(z);
alpha = alpha(idx);  
sidenum = sidenum(idx);
N = length(z);

% === Construction of target region ===

% Insert "anchor" points after every infinite vertex.
% (More than is really needed, but helpful for various boring reasons.)
add = find( alpha(1:N-1)==0 );
dz = min(z(add+1),z(N-1)+2) - z(add);
new = z(add) + dz/3;
new = [new; z(add) + 2*dz/3];
% Don't let the first vertex be infinite.
if alpha(N)==0
  new = [z(1)-2;z(1)-1;new]; 
end
z = [ z; new ];
alpha = [ alpha; ones(size(new)) ];
sidenum = [ sidenum; NaN*ones(size(new)) ];

% Sort again.
[z,idx] = sort(z);
alpha = alpha(idx);  sidenum = sidenum(idx);
N = length(z);

% Find side numbers.
original = find(~isnan(sidenum));
sidenum(1:original(1)-1) = n;
for j = 1:n-1
  sidenum( original(j)+1:original(j+1)-1 ) = j;
end

% Now that we have prevertices, transform to the disk since it's more
% convenient overall.
phi = moebius(z(N-2:N),[-1 -i 1]);
z = phi(z); 
zb = phi(zb);
beta = alpha-1;
qdata = scqdata(beta,15);

% We need to determine the multiplicative constant, using a known
% difference in real parts (D values). 
% First, find all finite D vertices.
Dvert = find(dirichlet(sidenum) & alpha~=0);
first = Dvert(1);
second = Dvert( 1 + min( find( diff( bdata(sidenum(Dvert)) ) ) ) ); 
% Use multiple versions and overdetermine for safety.
dt = angle( z(rem(first,N)+1)/z(first) );
z1 = z(first) * exp(i*dt*(0.2:0.2:0.8)');
dt = angle( z(rem(second,N)+1)/z(second) );
z2 = z(second) * exp(i*dt*(0.2:0.2:0.8)');
q = rsquad(z1,0*z1,0*z1,z,beta,zb,qdata) - ...
    rsquad(z2,0*z2,0*z1,z,beta,zb,qdata);
dphi = diff( bdata( sidenum([first second]) ) );
C = [ real(q) -imag(q) ] \ (dphi*ones(4,1));
c = C(1) + 1i*C(2);

% Start the target with known infinities.
w = NaN*ones(N,1);
w(alpha==0) = Inf;
finite = find( alpha > 0 );

% Integrate to find the finite ones, relative to the first.
idxl = finite(1:end-1)';  idxr = finite(2:end)';
zl = z(idxl); zr = z(idxr);
midpt = 0*zl;
q1 = rsquad(zl,midpt,idxl,z,beta,zb,qdata);
q2 = rsquad(zr,midpt,idxr,z,beta,zb,qdata);
w(finite) = c*cumsum([0;q1-q2]);

% Translate to correct first Dirichlet value. 
w = w - w(first) + bdata(sidenum(first));
Q = polygon(w,alpha);

% Find branch points.
z1 = repmat(z(1),[length(zb) 1]); 
branch = w(1) + c*rsquad(z1,zb,ones(size(z1)),z,beta,zb,qdata);

% Create output objects.
fd = diskmap(p,z(original));
if isempty(zb)
  fr = diskmap(Q,z,c);
else
  fr = riesurfmap(Q,branch,z,zb,c);
end
phi = composite( inv(fd), fr, inline('real(z)') );


% === Needed for RIESURFMAP quadrature. ===
function I = rsquad(z1,z2,varargin) 
%RSQUAD  (not intended for calling directly by the user)
%   Numerical quadrature for the Riemann surface map.

%   Copyright 2002 by Toby Driscoll.
%   $Id: lapsolve.m 298 2009-09-15 14:36:37Z driscoll $

%   RSQUAD(z1,z2,sing1,z,beta,zb,qdat)
%
%   z1,z2 are vectors of left and right endpoints.  sing1 is a vector of
%   integer indices which label the singularities in z1.  So if sing1(5)
%   = 3, then z1(5) = z(3).  A zero means no singularity.  z is the
%   vector of singularities; beta is the vector of associated turning
%   angles.  qdat is quadrature data from SCQDATA.
%
%   Make sure that z and beta are column vectors.
% 
%   The integral is subdivided, if necessary, so that no singularity
%   lies closer to the left endpoint than 1/2 the length of the
%   integration (sub)interval.
%
%   RSQUAD(z1,z2,sing1,sing2,z,beta,zb,qdat)
%  
%   Integrate from one singularity to another. Picks the origin as the
%   midpoint of integration; this may NOT be wise if the path goes near a
%   branch point.

if nargin==8
  % Break into two pieces with recursive call.
  [sing1,sing2,z,beta,zb,qdat] = deal(varargin{:});
  mid = zeros(size(z1));
  I1 = rsquad(z1,mid,sing1,z,beta,qdat);
  I2 = rsquad(z2,mid,sing2,z,beta,qdat);
  I = I1-I2;
  return
else
  [sing1,z,beta,zb,qdat] = deal(varargin{:});
end

nqpts = size(qdat,1);
n = length(z);
bigz = z(:,ones(1,nqpts));
bigbeta = beta(:,ones(1,nqpts));
if isempty(sing1)
  sing1 = zeros(length(z1),1);
end
B = length(zb);
if B>0
  bigbranch = zb(:,ones(1,nqpts));
else
  bigbranch = zeros(1,nqpts);
end

I = zeros(size(z1));
nontriv = find(z1(:)~=z2(:))';

for k = nontriv
  za = z1(k);
  zb = z2(k);
  sng = sing1(k);

  % Allowable integration step, based on nearest singularity.
  dist = min(1,2*min(abs(z([1:sng-1,sng+1:n])-za))/abs(zb-za));
  zr = za + dist*(zb-za);
  % Adjust Gauss-Jacobi nodes and weights to interval.
  ind = rem(sng+n,n+1)+1;
  nd = ((zr-za)*qdat(:,ind) + zr + za).'/2; % G-J nodes
  wt = ((zr-za)/2) * qdat(:,ind+n+1);	% G-J weights
  terms = 1 - (nd(ones(n,1),:))./bigz;
  if any(terms(:)==0)
    % Endpoints are practically coincident.
    I(k) = 0;
  else
    % Use Gauss-Jacobi on first subinterval, if necessary.
    if sng > 0
      terms(sng,:) = terms(sng,:)./abs(terms(sng,:));
      wt = wt*(abs(zr-za)/2)^beta(sng);
    end
    Q = exp(sum(log(terms).*bigbeta));
    if B > 0
      ND = nd(ones(B,1),:);
      Q = Q.*prod((ND-bigbranch).*(1-ND.*conj(bigbranch)),1);
    end
    I(k) = Q*wt;
    while dist < 1              
      % Do regular Gaussian quad on other subintervals.
      zl = zr;
      dist = min(1,2*min(abs(z-zl))/abs(zl-zb));
      zr = zl + dist*(zb-zl);
      nd = ((zr-zl)*qdat(:,n+1) + zr + zl).'/2;
      wt = ((zr-zl)/2) * qdat(:,2*n+2);
      Q = exp(sum(log(1 - nd(ones(n,1),:)./bigz).*bigbeta));
      if B > 0
        ND = nd(ones(B,1),:);
        Q = Q.*prod((ND-bigbranch).*(1-ND.*conj(bigbranch)),1);
      end      
      I(k) = I(k) + Q*wt;
    end
  end
end

