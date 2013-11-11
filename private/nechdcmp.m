function [L,maxadd] = nechdcmp(H,maxoffl)
%
% [L,maxadd] = nechdcmp(H,maxoffl)
%
% This function is part of the Nonlinear Equations package, see NESOLVE.M.
%
% This is a "Perturbed Cholesky Decomposition".  It finds a lower
% triangular matrix L such that LL' is a factorization of H+D, where
% D is a diagonal (non-negative) matrix that is added to H if necessary
% to make it positive definite (so that the factorization is possible).
% If H is already positive definite, the ordinary Cholesky decomposition
% (D=0) is carried out.
%
% Algorithm A5.5.2: Part of the modular software system from
% the appendix of the book "Numerical Methods for Unconstrained
% Optimization and Nonlinear Equations" by Dennis & Schnabel 1983.
%
% Coded in Matlab by Sherkat Masoum M., April 1988.
% Edited by Richard T. Behrens, June 1988.
%

%
% Check input arguments.
%
[m,n]=size(H);
if (m ~=n)
  error('Matrix H must be square.')
end

%
% Algorithm step 1.
%
minl=(eps^.25) * maxoffl;

%
% Algorithm step 2.
%
if (maxoffl == 0.)
% This is the case when H is known to be positive def.
  maxoffl=sqrt(max(diag(H)));
  minl2=(eps^.5) * maxoffl;
end

%
% Algorithm step 3.
%
maxadd=0.;      % the maximum diagonal element (so far) in D.

%
% Algorithm step 4.
%
for j=1:n
  if (j==1)
    L(j,j)=H(j,j);
  else
    L(j,j)=H(j,j)-L(j,1:j-1)*L(j,1:j-1)';
  end
  minljj=0.;
  for i=j+1:n
    if (j==1)
      L(i,j)=H(j,i);
    else
      L(i,j)=H(j,i)-L(i,1:j-1)*L(j,1:j-1)';
    end
    minljj=max(abs(L(i,j)),minljj);
  end
  minljj=max(minljj/maxoffl,minl);
  if (L(j,j) > minljj^2)
    % Normal Cholesky iteration
    L(j,j)=sqrt(L(j,j));
  else
    % Augment H(j,j)
    if (minljj < minl2)
      minljj=minl2;
      % Only possible when input maxoffl=0
    end
    maxadd=max(maxadd,(minljj^2-L(j,j)));
    L(j,j)=minljj;
  end
  for i=j+1:n
    L(i,j)=L(i,j)/L(j,j);
  end
  %L(j+1:n,j)=L(j+1:n,j)/L(j,j);
end

