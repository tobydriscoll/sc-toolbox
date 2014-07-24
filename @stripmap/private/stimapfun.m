function zdot = stimapfun(wp,yp,flag,scale,z,beta,c);

%   Used by STINVMAP for solution of an ODE.

%   Copyright 1998 by Toby Driscoll.
%   $Id: stimapfun.m 298 2009-09-15 14:36:37Z driscoll $

lenyp = length(yp);
lenzp = lenyp/2;
zp = yp(1:lenzp) + i*yp(lenzp+1:lenyp);

f = scale./stderiv(zp,z,beta,c);
zdot = [real(f);imag(f)];
