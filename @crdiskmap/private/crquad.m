function I = crquad(z1,sing1,z,beta,qdat)
%CRQUAD (not intended for calling directly by the user)
%   Numerical quadrature for the disk map, crossratio formulation.

%   z1 is a vector of left endpoints; 0 is always the right endpoint.
%   sing1 is a vector of integer indices which label the singularities
%   in z1.  So if sing1(5) = 3, then z1(5) = z(3).  A zero means no
%   singularity.  z is the vector of prevertices; beta is the vector of
%   associated turning angles.  qdat is quadrature data from SCQDATA.
%
%   The integral is subdivided, if necessary, so that no singularity
%   lies closer to the left endpoint than 1/2 the length of the
%   integration (sub)interval.
%
%   The conceptual differences between CRQUAD and DQUAD are tiny.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crquad.m 212 2002-09-25 17:31:37Z driscoll $

nqpts = size(qdat,1);
n = length(z);
bigz = z(:,ones(1,nqpts));
bigbeta = beta(:,ones(1,nqpts));
if isempty(sing1)
  sing1 = zeros(length(z1),1);
end

% Will ignore zero values of beta
ignore = abs(beta) < eps;

I = zeros(size(z1));
nontriv = find(abs(z1(:))>eps)';

for k = nontriv
  za = z1(k);
  sng = sing1(k);
  mask = ~ignore;

  % Integration subintervals are based on nearest singularity
  if sng
    mask(sng) = 0;
  end
  mindist = min(abs(z(mask)-za));
  panels = max(1,ceil(-log(mindist)/log(2)));
  if sng
    mask(sng) = ~ignore(sng);
  end

  qcol = rem(sng+n,n+1)+1;
  keep = find(mask);

  % Adjust sng for ignored vertices
  if sng
    if ignore(sng)
      sng = 0;
    else
      sng = sng - sum(ignore(1:sng));
    end
  end
  
  for j=1:panels
    if j==1
      h = 2^(1-panels);
    else
      h = h + 2^(j-panels-1);
      za = zb;
      qcol = n+1;
    end    
    zb = z1(k)*(1-h);
  
    % Adjust Gauss-Jacobi nodes and weights to interval.
    nd = ((zb-za)*qdat(:,qcol) + zb + za)/2; % G-J nodes
    wt = ((zb-za)/2) * qdat(:,qcol+n+1);     % G-J weights
    terms = 1 - (nd(:,ones(sum(mask),1)).')./bigz(mask,:);
    % Check for coincident values indicating crowding. (Should never
    % happen!)
    if any( diff(nd)==0 ) || any(terms(:)==0)
      warning('Prevertices are too crowded.')
      I(k) = 0;
    else
      % Use Gauss-Jacobi on first subinterval, if necessary.
      if sng && (qcol < n+1)
	terms(sng,:) = terms(sng,:)./abs(terms(sng,:));
	wt = wt*(abs(zb-za)/2)^beta(keep(sng));
      end
      I(k) = I(k) + exp(sum(log(terms).*bigbeta(mask,:)))*wt;
    end
  end
end

