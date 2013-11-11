function wp = eval(M,zp,tol)
%EVAL Evaluate Schwarz-Christoffel crossratio disk map at points.
%   EVAL(M,ZP) evaluates the Schwarz-Christoffel map M at the points ZP
%   in the unit disk. The default tolerance of M is used.
%   
%   EVAL(M,ZP,TOL) attempts to give an answer accurate to TOL. If TOL is
%   less than the accuracy of M, this is unlikely to be met.
%   
%   See also CRDISKMAP, EVALINV.

%   Copyright 1998 by Toby Driscoll.
%   $Id: eval.m 7 1998-05-10 04:37:19Z tad $

if nargin < 3
  qdata = M.qdata;
else 
  qdata = tol;
end

p = polygon(M);
w = vertex(p);
beta = angle(p) - 1;
cr = M.crossratio;
aff = M.affine;
wcfix = M.center{2};
Q = M.qlgraph;

wp = NaN*zp;
idx = abs(zp) <= 1+eps;
wp(idx) = crmap(zp(idx),w,beta,cr,aff,wcfix,Q,qdata);
