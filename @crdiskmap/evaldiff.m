function fp = evaldiff(M,zp)
%EVALDIFF Derivative of Schwarz-Christoffel crossratio disk map at points.
%   EVALDIFF(M,ZP) computes the derivative of the Schwarz-Christoffel
%   disk map M at the points ZP.
%   
%   See also CRDISKMAP, EVAL.

%   Copyright 1998 by Toby Driscoll.
%   $Id: evaldiff.m 7 1998-05-10 04:37:19Z tad $

beta = angle(polygon(M)) - 1;
cr = M.crossratio;
aff = M.affine;
wcfix = M.center{2};
Q = M.qlgraph;
qdata = M.qdata;

fp = crderiv(zp,beta,cr,aff,wcfix,Q,qdata);
