function fp = evaldiff(f,zp)
%EVALDIFF Evaluate derivative of Schwarz-Christoffel Riemann surface map.
%   EVALDIFF(F,ZP) computes the derivative of the Schwarz-Christoffel
%   Riemann surface map F at the points ZP.
%   
%   See also RIESURFMAP, EVAL.

%   Copyright 2002 by Toby Driscoll.
%   $Id: evaldiff.m 203 2002-09-13 20:09:37Z driscoll $

z = f.prevertex;
c = f.constant;
beta = angle(polygon(f)) - 1;
zb = f.prebranch;

fp = rsderiv(zp,z,beta,zb,c);
