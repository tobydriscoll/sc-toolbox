function wp = eval(M,zp,tol)
%EVAL Evaluate Schwarz-Christoffel half-plane map at points.
%   EVAL(M,ZP) evaluates the Schwarz-Christoffel map M at the points
%   ZP in the upper half-plane. The default tolerance of M is used.
%   
%   EVAL(M,ZP,TOL) attempts to give an answer accurate to TOL. If TOL
%   is less than the accuracy of M, this is unlikely to be met.
%   
%   See also HPLMAP, EVALINV.

%   Copyright 1998 by Toby Driscoll.
%   $Id: eval.m 7 1998-05-10 04:37:19Z tad $

p = polygon(M);
n = length(p);

if nargin < 3
  qdata = M.qdata;
else 
  qdata = tol;
end

wp = NaN*zp;
idx = imag(zp) > -eps;
wp(idx) = hpmap(zp(idx),vertex(p),angle(p)-1,M.prevertex,M.constant,qdata);
