function I = rsquad(z1,z2,varargin) 
%RSQUAD  (not intended for calling directly by the user)
%   Numerical quadrature for the Riemann surface map.

%   Copyright 2002 by Toby Driscoll.
%   $Id: rsquad.m 298 2009-09-15 14:36:37Z driscoll $

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
  I1 = rsquad(z1,mid,sing1,z,beta,zb,qdat);
  I2 = rsquad(z2,mid,sing2,z,beta,zb,qdat);
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

