function wp = eval(M,zp,tol)
%EVAL Evaluate Schwarz-Christoffel exterior map at points.
%   EVAL(M,ZP) evaluates the Schwarz-Christoffel map M at the points
%   ZP in the unit disk. The default tolerance of M is used.
%   
%   EVAL(M,ZP,TOL) attempts to give an answer accurate to TOL. If TOL
%   is less than the accuracy of M, this is unlikely to be met.
%   
%   See also EXTERMAP, EVALINV.

%   Copyright 1998 by Toby Driscoll.
%   $Id: eval.m 7 1998-05-10 04:37:19Z tad $

if nargin < 3
  qdata = M.qdata;
else 
  qdata = tol;
end

p = polygon(M);
w = flipud(vertex(p));
beta = flipud(1 - angle(p));
n = length(w);

wp = NaN*zp;
idx = abs(zp) <= 1+eps;
wp(idx) = demap(zp(idx),w,beta,M.prevertex,M.constant,qdata);
