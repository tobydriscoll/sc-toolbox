function fp = evaldiff(M,zp)
%EVALDIFF Derivative of Schwarz-Christoffel exterior map at points.
%   EVALDIFF(M,ZP) computes the derivative of the Schwarz-Christoffel
%   exterior map M at the points ZP.
%   
%   See also EXTERMAP, EVAL.

%   Copyright 1998 by Toby Driscoll.
%   $Id: evaldiff.m 7 1998-05-10 04:37:19Z tad $

z = M.prevertex;
c = M.constant;
beta = flipud(1 - angle(polygon(M)));

fp = dederiv(zp,z,beta,c);
