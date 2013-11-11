function I = dequad(z1,z2,sing1,z,beta,qdat)
%DEQUAD (not intended for calling directly by the user)
%   Numerical quadrature for the exterior map.

%   Copyright 1998 by Toby Driscoll.
%   $Id: dequad.m 212 2002-09-25 17:31:37Z driscoll $

%   z1,z2 are vectors of left and right endpoints.  sing1 is a vector of
%   integer indices which label the singularities in z1.  So if sing1(5)
%   = 3, then z1(5) = z(3).  A zero means no singularity.  z is the
%   vector of prevertices (all singularities except the origin); beta is
%   the vector of associated turning angles.  qdat is quadrature data
%   from SCQDATA.
%
%   Make sure that z and beta are column vectors.
%   
%   DEQUAD integrates from a possible singularity at the left end to a
%   regular point at the right.  If both endpoints are singularities,
%   you must break the integral into two pieces and make two calls.
%   
%   The integral is subdivided, if necessary, so that no singularity
%   lies closer to the left endpoint than 1/2 the length of the
%   integration (sub)interval.  But the singularity at the origin is NOT
%   accounted for in this decision.

nqpts = size(qdat,1);
n = length(z);
bigz = z(:,ones(1,nqpts));
beta = [beta(:);-2];
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
  ind = sng + (n+1)*(sng==0);
  nd = ((zr-za)*qdat(:,ind) + zr + za).'/2; % nodes
  wt = ((zr-za)/2) * qdat(:,ind+n+1);	% weights
  terms = 1 - nd(ones(n,1),:)./bigz;
  if any(terms(:)==0)
    % Endpoints are practically coincident.
    I(k) = 0;
  else 
    terms = [terms;nd];
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
      terms = 1 - nd(ones(n,1),:)./bigz;
      terms = [terms;nd];
      I(k) = I(k) + exp(sum(log(terms).*bigbeta))*wt;
    end
  end
end

