function I = dquad(z1,z2,sing1,z,beta,qdat)
%DQUAD: Numerical quadrature for the disk map.
%   DQUAD(ZL,ZR,S,Z,BETA,QDAT,MIDPT) performs the integration for the SC
%   disk formula. ZL and ZR are vectors of left and right endpoints. S
%   is a vector of integers. If ZL(k) = Z(m), then S(k) should have the
%   value m; otherwise, S(k) should be zero. 
%   
%   Z and BETA are prevertices and turning angles for the SC map. QDAT
%   is a matrix of quadrature data (see SCQDATA). 
%   
%   The integration is adaptive in the sense that members of Z (with
%   nonzero BETA) that are close to the left endpoints cause
%   subdivision. This is NOT true of singularities close to the right end.

%   Copyright 1998--2001 by Toby Driscoll.
%   $Id: dquad.m 212 2002-09-25 17:31:37Z driscoll $

%   z1,z2 are vectors of left and right endpoints.  sing1 is a vector of
%   integer indices which label the singularities in z1.  So if sing1(5)
%   = 3, then z1(5) = z(3).  A zero means no singularity.  z is the
%   vector of singularities; beta is the vector of associated turning
%   angles.  qdat is quadrature data from SCQDATA.
%
%   Make sure that z and beta are column vectors.
%   
%   DQUAD integrates from a possible singularity at the left end to a
%   regular point at the right.  If both endpoints are singularities,
%   you must break the integral into two pieces and make two calls.
%   
%   The integral is subdivided, if necessary, so that no singularity
%   lies closer to the left endpoint than 1/2 the length of the
%   integration (sub)interval.

nqpts = size(qdat,1);
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
    I(k) = exp(sum(log(terms).*bigbeta))*wt;
    while dist < 1              
      % Do regular Gaussian quad on other subintervals.
      zl = zr;
      dist = min(1,2*min(abs(z-zl))/abs(zl-zb));
      zr = zl + dist*(zb-zl);
      nd = ((zr-zl)*qdat(:,n+1) + zr + zl).'/2;
      wt = ((zr-zl)/2) * qdat(:,2*n+2);
      I(k) = I(k) + exp(sum(log(1 - nd(ones(n,1),:)./bigz).*bigbeta)) * wt;
    end
  end
end
