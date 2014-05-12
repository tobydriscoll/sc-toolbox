function I = hpquad(z1,z2,varargin)
%HPQUAD (not intended for calling directly by the user)
%   Numerical quadrature for the half-plane map.

%   Copyright 1998 by Toby Driscoll.
%   $Id: hpquad.m 298 2009-09-15 14:36:37Z driscoll $

%   HPQUAD(z1,z2,sing1,z,beta,qdat)
%   z1,z2 are vectors of left and right endpoints.  sing1 is a vector of
%   integer indices which label the singularities in z1.  So if sing1(5)
%   = 3, then z1(5) = z(3).  A zero means no singularity.  z is the
%   vector of finite singularities; beta is the vector of associated
%   turning angles.  qdat is quadrature data from SCQDATA.
%
%   Make sure z and beta are column vectors.
%   
%   HPQUAD integrates from a possible singularity at the left end to a
%   regular point at the right.  If both endpoints are singularities,
%   you must break the integral into two pieces and make two calls, or
%   call HPQUAD(z1,z2,sing1,sing2,z,beta,qdat) and accept an automatic
%   choice. 
%   
%   The integral is subdivided, if necessary, so that no singularity
%   lies closer to the left endpoint than 1/2 the length of the
%   integration (sub)interval.

if nargin==7
  % Break into two pieces with recursive call.
  [sing1,sing2,z,beta,qdat] = deal(varargin{:});
  mid = (z1+z2)/2;
  mid = mid + 1i*abs(mid);
  I1 = hpquad(z1,mid,sing1,z,beta,qdat);
  I2 = hpquad(z2,mid,sing2,z,beta,qdat);
  I = I1-I2;
  return
else
  [sing1,z,beta,qdat] = deal(varargin{:});
end

nqpts = size(qdat,1);
% Note: Here n is the total number of *finite* singularities; i.e., the
% number of terms in the product appearing in the integrand.
n = length(z);
bigz = z(:,ones(1,nqpts));
bigbeta = beta(:,ones(1,nqpts));
if isempty(sing1)
  sing1 = zeros(length(z1),1);
end

I = zeros(size(z1));
nontriv = find(z1(:)~=z2(:))';

for k = nontriv
  za = z1(k);
  zb = z2(k);
  sng = sing1(k);

  % Allowable integration step, based on nearest singularity.
  dist = min(1,2*min(abs(z([1:sng-1,sng+1:n])-za))/abs(zb-za));
%%  if isempty(dist), dist=1; end
  zr = za + dist*(zb-za);
  ind = rem(sng+n,n+1)+1;
  % Adjust Gauss-Jacobi nodes and weights to interval.
  nd = ((zr-za)*qdat(:,ind) + zr + za).'/2; % G-J nodes
  wt = ((zr-za)/2) * qdat(:,ind+n+1); 	% G-J weights
  terms = nd(ones(n,1),:) - bigz;
  if any(terms(:)==0) 
    % Endpoints are practically coincident.
    I(k) = 0;
  else
    % Use Gauss-Jacobi on first subinterval, if necessary.
    if sng > 0
      terms(sng,:) = terms(sng,:)./abs(terms(sng,:));
      wt = wt*(abs(zr-za)/2)^beta(sng);
    end
    I(k) = exp(sum(log(terms).*bigbeta,1))*wt;
    while dist < 1              
      % Do regular Gaussian quad on other subintervals.
      zl = zr;
      dist = min(1,2*min(abs(z-zl))/abs(zl-zb));
      zr = zl + dist*(zb-zl);
      nd = ((zr-zl)*qdat(:,n+1) + zr + zl).'/2;
      wt = ((zr-zl)/2) * qdat(:,2*n+2);
      terms = nd(ones(n,1),:) - bigz;
      I(k) = I(k) + exp(sum(log(terms).*bigbeta,1)) * wt;
    end
  end
end

