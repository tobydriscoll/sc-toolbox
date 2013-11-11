function f = eval(map,z)
%Evaluate Moebius transformation at point(s).
%   EVAL(M,Z) evaluates the Moebius transformation M at the point(s)
%   in Z. Infinity is a valid input.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: eval.m 93 2000-05-24 17:57:50Z tad $

f = NaN*zeros(size(z));
atinf = isinf(z);
if any(atinf)
  f(atinf) = map.coeff(2)/map.coeff(4);
end

num = map.coeff(2)*z(~atinf) + map.coeff(1);
den = map.coeff(4)*z(~atinf) + map.coeff(3);

toinf = abs(den) < 3*eps;
den(toinf) = NaN;
num(toinf) = 1;

f(~atinf) = num./den;
f(isnan(f)) = Inf;