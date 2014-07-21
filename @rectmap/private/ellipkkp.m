function [K,Kp] = ellipkkp(L)
%ELLIPKKP Complete elliptic integral of the first kind, with complement.
%   K = ELLIPKKP(L) returns the value of the complete elliptic
%   integral of the first kind, evaluated at M=exp(-2*pi*L), 0 < L
%   < Inf.
%   
%   [K,KP] = ELLIPKKP(L) also returns the result for complementary
%   parameter 1-M, which is useful when M < EPS.  Even when M <
%   1e-6, the built-in ELLIPKE can lose digits of accuracy for KP.
% 
%   Recall that the elliptic modulus k is related to the parameter
%   M by M = k^2.
% 
%   Copyright (c)1999 by Toby Driscoll. 
%   $Id: ellipkkp.m 298 2009-09-15 14:36:37Z driscoll $ 

%   ELLIPKKP uses the method of the arithmetic-geometric mean described
%   in 17.6 of M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%   Functions," Dover, 1965.  Same method as in ELLIPKE, only
%   interchanging 1 and 1-m to find KP.

% When m=exp(-2*pi*L) is extremely small, use O(m) approximations.
if L > 10
  K = pi/2;
  Kp = pi*L + log(4);
  return
end

m = exp(-2*pi*L);
a0 = 1;
b0 = sqrt(1-m);
s0 = m;
i1 = 0; mm = 1;
while mm > eps
  a1 = (a0+b0)/2;
  b1 = sqrt(a0.*b0);
  c1 = (a0-b0)/2;
  i1 = i1 + 1;
  w1 = 2^i1*c1.^2;
  mm = max(max(w1));
  s0 = s0 + w1;
  a0 = a1;
  b0 = b1;
end
K = pi./(2*a1);

im = find(m==1);
if ~isempty(im)
  K(im) = K(im)*inf;
end

if nargout > 1
  a0 = 1;
  b0 = sqrt(m);
  s0 = 1-m;
  i1 = 0; mm = 1;
  while mm > eps
    a1 = (a0+b0)/2;
    b1 = sqrt(a0.*b0);
    c1 = (a0-b0)/2;
    i1 = i1 + 1;
    w1 = 2^i1*c1.^2;
    mm = max(max(w1));
    s0 = s0 + w1;
    a0 = a1;
    b0 = b1;
  end
  Kp = pi./(2*a1);
  im = find(m==0);
  if ~isempty(im)
    Kp(im) = Kp(im)*inf;
  end
end

