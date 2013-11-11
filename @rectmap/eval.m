function wp = eval(M,zp,tol)
%EVAL Evaluate Schwarz-Christoffel rectangle map at points.
%   EVAL(M,ZP) evaluates the Schwarz-Christoffel map M at the points ZP
%   in the source rectangle of M. The default tolerance of M is used.
%   
%   EVAL(M,ZP,TOL) attempts to give an answer accurate to TOL. If TOL is
%   less than the accuracy of M, this is unlikely to be met.
%   
%   See also RECTMAP, EVALINV.

%   Copyright 1998 by Toby Driscoll.
%   $Id: eval.m 267 2003-05-08 18:07:51Z driscoll $

p = polygon(M);
n = length(p);
z = M.prevertex;
zr = z(corners(M));

if nargin < 3
  qdata = M.qdata;
else 
  qdata = tol;
end

wp = NaN*zp;
opt = scmapopt(M);
idx = find( isinpoly(zp,polygon(zr),opt.Tolerance) );
wp(idx) = ...
    rmap(zp(idx),vertex(p),angle(p)-1,z,M.constant,M.stripL,qdata);
