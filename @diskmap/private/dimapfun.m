function zdot = dimapfun(wp,yp,scale,z,beta,c)
%   Used by DINVMAP for solution of an ODE.

%   Copyright 1998 by Toby Driscoll.
%   $Id: dimapfun.m 7 1998-05-10 04:37:19Z tad $

lenyp = length(yp);
lenzp = lenyp/2;
zp = (yp(1:lenzp)+sqrt(-1)*yp(lenzp+1:lenyp));

f = scale./dderiv(zp,z,beta,c);
zdot = [real(f);imag(f)];

end
