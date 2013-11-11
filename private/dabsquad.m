function I = dabsquad(z1,z2,sing1,z,beta,qdat)
%DABSQUAD (not intended for calling directly by the user)
%   Numerical quadrature for side lengths in the disk map.
%
%   DABSQUAD is similar to DQUAD, but the integrand is replaced by its
%   absolute value and the path of integration is along the unit circle
%   rather than the straight line between endpoints. This allows the use
%   of real rather than complex logarithms in the computation of powers,
%   which can result in substantial savings.

%   z1,z2 are vectors of left and right endpoints.  sing1 is a vector of
%   integer indices which label the singularities in z1.  So if sing1(5)
%   = 3, then z1(5) = z(3).  A zero means no singularity.  z is the
%   vector of singularities; beta is the vector of associated turning
%   angles.  qdat is quadrature data from SCQDATA.
%
%   Make sure that z and beta are column vectors.
%   
%   DABSQUAD integrates from a possible singularity at the left end to a
%   regular point at the right.  If both endpoints are singularities,
%   you must break the integral into two pieces and make two calls. The
%   integration always takes along the smaller arc of the circle between
%   the endpoints.
%   
%   The integral is subdivided, if necessary, so that no singularity
%   lies closer to the left endpoint than 1/2 the length of the
%   integration (sub)interval.

%   Copyright 1997 by Toby Driscoll.  Last updated 05/08/97.

nqpts = size(qdat,1);
n = length(z);
argz = angle(z);
argz1 = angle(z1);
argz2 = angle(z2);
ang21 = angle(z2./z1);
% Modification needed if integration passes through -1
discont = (argz2-argz1) .* ang21 < 0;
argz2(discont) = argz2(discont) + 2*pi*sign(ang21(discont));
bigargz = argz(:,ones(1,nqpts));
bigbeta = beta(:,ones(1,nqpts));
if isempty(sing1)
  sing1 = zeros(length(z1),1);
end

I = zeros(size(z1));
nontriv = find(z1(:)~=z2(:))';

for k = nontriv
  arga = argz1(k);
  argb = argz2(k);
  za = z1(k);
  zb = z2(k);
  sng = sing1(k);

  % Allowable integration step, based on nearest singularity.
  dist = min(1,2*min(abs(z([1:sng-1,sng+1:n])-za))/abs(zb-za));
  argr = arga + dist*(argb-arga);
  % Adjust Gauss-Jacobi nodes and weights to interval.
  ind = rem(sng+n,n+1)+1;
  nd = ((argr-arga)*qdat(:,ind) + argr + arga).'/2; % G-J nodes
  wt = (abs(argr-arga)/2)*qdat(:,ind+n+1);	% G-J weights
  theta = rem(nd(ones(n,1),:)-bigargz+2*pi, 2*pi);
  mask = theta > pi;
  theta(mask) = 2*pi - theta(mask);
  terms = 2*sin(theta/2);
  if any(any(~terms))
    % Endpoints are practically coincident.
    I(k) = 0;
  else
    % Use Gauss-Jacobi on first subinterval, if necessary.
    if sng > 0
      terms(sng,:) = terms(sng,:)./abs(nd - arga);
      wt = wt*(abs(argr-arga)/2)^beta(sng);
    end
    I(k) = exp(sum(log(terms).*bigbeta))*wt;
    while dist < 1              
      % Do regular Gaussian quad on other subintervals.
      argl = argr;
      zl = exp(i*argl);
      dist = min(1,2*min(abs(z-zl)/abs(zl-zb)));
      argr = argl + dist*(argb-argl);
      nd = ((argr-argl)*qdat(:,n+1) + argr + argl).'/2;
      wt = (abs(argr-argl)/2) * qdat(:,2*n+2);
      theta = rem(nd(ones(n,1),:)-bigargz+2*pi, 2*pi);
      theta(theta>pi) = 2*pi - theta(theta>pi);
      I(k) = I(k) + exp(sum(log(2*sin(theta/2)).*bigbeta)) * wt;
    end
  end
end
