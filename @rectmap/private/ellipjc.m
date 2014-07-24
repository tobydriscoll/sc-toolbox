function [sn,cn,dn] = ellipjc(u,L,flag)
%ELLIPJC Jacobi elliptic functions for complex argument.
%   [SN,CN,DN] = ELLIPJC(U,L) returns the values of the Jacobi
%   elliptic functions evaluated at complex argument U and
%   parameter M=exp(-2*pi*L), 0 < L < Inf.  Recall that M = k^2,
%   where k is the elliptic modulus.
%   
%   U may be a matrix; L must be a scalar.  The entries of U are
%   expected to lie within the rectangle |Re U| < K, 0 < Im U <
%   Kp, where [K,Kp] = ELLIPK(L).
%   
%   Copyright (c) 1999 by Toby Driscoll. 
%   $Id: ellipjc.m 298 2009-09-15 14:36:37Z driscoll $

%   The built-in ELLIPJ can't handle compelx arguments, and
%   standard transformations to handle this would require ELLIPJ
%   called with parameter 1-M. When M < eps (or is even close),
%   this can't be done accurately.
%   
%   The algorithm is the descending Landen transformation,
%   described in L. Howell's PhD thesis from MIT. Additional
%   formulas from Gradshteyn & Ryzhik, 5th ed., and Abramowitz
%   & Stegun.

if nargin < 3
  % Absence of flag parameter indicates we must check for and transform u in
  % the upper half of the rectangle.
  [K,Kp] = ellipkkp(L);
  high = imag(u) > Kp/2;
  u(high) = i*Kp - u(high);
  m = exp(-2*pi*L);
else
  % Recursive call--L is actually m.
  high = zeros(size(u));
  m = L;
end

if m < 4*eps
  sinu = sin(u);
  cosu = cos(u);
  sn = sinu + m/4*(sinu.*cosu-u).*cosu;
  cn = cosu + m/4*(-sinu.*cosu+u).*sinu;
  dn = 1 + m/4*(cosu.^2-sinu.^2-1);
else
  if m > 1e-3
    kappa = (1-sqrt(1-m))/(1+sqrt(1-m));
  else
    kappa = polyval([132,42,14,5,2,1,0],m/4);
  end
  mu = kappa^2;
  v = u/(1+kappa);
  [sn1,cn1,dn1] = ellipjc(v,mu,1);
  denom = (1+kappa*sn1.^2);
  sn = (1+kappa)*sn1 ./ denom;
  cn = cn1.*dn1 ./ denom;
  dn = (1-kappa*sn1.^2) ./ denom;
end

if any(high(:))
  snh = sn(high);
  cnh = cn(high);
  dnh = dn(high);
  sn(high) = -1./(sqrt(m)*snh);
  cn(high) = i*dnh./(sqrt(m)*snh);
  dn(high) = i*cnh./snh;
end
