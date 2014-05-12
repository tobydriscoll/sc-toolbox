function zdot = rimapfun(wp,yp,flag,scale,z,beta,c,zs,L);
%   Used by RINVMAP for solution of an ODE.

%   Copyright 1998 by Toby Driscoll.
%   $Id: rimapfun.m 298 2009-09-15 14:36:37Z driscoll $

lenyp = length(yp);
lenzp = lenyp/2;
zp = yp(1:lenzp) + i*yp(lenzp+1:lenyp);

f = scale./rderiv(zp,z,beta,c,L,zs);
zdot = [real(f);imag(f)];
