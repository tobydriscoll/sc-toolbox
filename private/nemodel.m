function [m,h,sn] = nemodel(fc,J,g,sf,sx,globmeth)
%
% [m,h,sn] = nemodel(fc,J,g,sf,sx,globmeth)
%
% This function is part of the Nonlinear Equations package, see NESOLVE.M.
%
% It forms the affine model for use in solving nonlinear equations.
%
% Algorithm A6.5.1: Part of the modular software system from
% the appendix of the book "Numerical Methods for Unconstrained
% Optimization and Nonlinear Equations" by Dennis & Schnabel 1983.
%
% Coded in Matlab by Sherkat Masoum M., March 1988.
% Edited by Richard T. Behrens, June 1988.
%
%
% Algorithm step 1.
%
n = length(J);
m = diag(sf)*J;

%
% Algorithm step 2.
%
[m,m1,m2,sing]=neqrdcmp(m);

%
% Algorithm step 3.
%
if (sing == 0)
  for j=2:n
    m(1:(j-1),j)=m(1:(j-1),j)/sx(j);
  end
  m2 = m2./sx;
  est = neconest(m,m2);
else
  est=0.;
end

%
% Algorithm step 4.
%
if (sing ==1) | (est > 1./eps) | isnan(est)
  h = J'*diag(sf);
  h = h*h';
  % calculate hnorm=norm(invDxHinvDx)
  tem = abs(h(1,:)) * (ones(n,1)./sx);
  hnorm=(1./sx(1))*tem;
  for i=2:n
    tem1=sum(abs(h(:,i))./sx);
    tem2=sum(abs(h(i,:))./(sx.'));
    temp=(1./sx(i))/(tem1+tem2);
    hnorm=max(temp,hnorm);
  end
  h = h + sqrt(n*eps) * hnorm * (diag(sx)^2);
  % caculate sn=inv(H)*g, and keep m (the cholesky factor) for later use.
  [m,maxadd] = nechdcmp(h,0);
  sn = -m'\(m\g);
else
  % Calculate normal Newton step
  for j=2:n
    m(1:(j-1),j)=m(1:(j-1),j)*sx(j);
  end
  m2 = m2.*sx;
  sn = -sf.*fc;
  sn = neqrsolv(m,m1,m2,sn);
  if (globmeth ==2) | (globmeth ==3)
    % the cholesky factor (for later use) is the same as R' from QR.
    m = triu(m) + triu(m)';
    m = m - diag(diag(m)) + diag(m2);
  end
  if (globmeth == 2)
    L = tril(m);
    h = L*L';          % This is J'*J, an approximation of H.
  else
    h = [];
  end
end

