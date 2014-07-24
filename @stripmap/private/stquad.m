function I = stquad(z1,z2,sing1,z,beta,qdat)
%STQUAD (not intended for calling directly by the user)
%   Numerical quadrature for the strip map.

%   Copyright 1998 by Toby Driscoll.
%   $Id: stquad.m 298 2009-09-15 14:36:37Z driscoll $

%   z1,z2 are vectors of left and right endpoints.  sing1 is a vector of
%   integer indices which label the singularities in z1.  So if sing1(5)
%   = 3, then z1(5) = z(3).  A zero means no singularity.  z is the
%   vector of *all* singularities, including the "ends" of the strip at
%   \pm Inf. beta is the vector of associated turning angles. qdat is
%   quadrature data from SCQDATA. It should include all the beta values,
%   even though the ends are never used in this manner.
%
%   Make sure z and beta are column vectors.
%   
%   STQUAD integrates from a possible singularity at the left end to a
%   regular point at the right.  If both endpoints are singularities,
%   you must break the integral into two pieces and make two calls.
%   
%   The integral is subdivided, if necessary, so that no singularity
%   lies closer to the left endpoint than 1/2 the length of the
%   integration (sub)interval.

n = length(z);
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
  ind = rem(sng+n,n+1)+1;
  % Adjust Gauss-Jacobi nodes and weights to interval.
  nd = ((zr-za)*qdat(:,ind) + zr + za).'/2; % G-J nodes
  wt = ((zr-za)/2) * qdat(:,ind+n+1); 	% G-J weights
  if any( diff([za;nd(:);zr])==0 ) 
    % Endpoints are practically coincident.
    I(k) = 0;
  else
    % Use Gauss-Jacobi on first subinterval, if necessary.
    if sng > 0
      wt = wt*(abs(zr-za)/2)^beta(sng);
    end
    I(k) = stderiv(nd,z,beta,1,sng)*wt;
    while (dist < 1) & ~isnan(I(k))
      % Do regular Gaussian quad on other subintervals.
      zl = zr;
      dist = min(1,2*min(abs(z-zl))/abs(zl-zb));
      zr = zl + dist*(zb-zl);
      nd = ((zr-zl)*qdat(:,n+1) + zr + zl).'/2;
      wt = ((zr-zl)/2) * qdat(:,2*n+2);
      I(k) = I(k) + stderiv(nd,z,beta,1)*wt;
    end
  end
end

