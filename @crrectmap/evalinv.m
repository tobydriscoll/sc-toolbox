function wp = evalinv(M,zp,tol)
%EVALINV Invert Schwarz-Christoffel crossratio rectified map at points.
%   EVALINV(M,ZP) evaluates the inverse of the Schwarz-Christoffel map M
%   at the points WP in the polygon. The default tolerance of M is used.
%   
%   EVALINV(M,ZP,TOL) attempts to give an answer accurate to TOL. If TOL
%   is less than the accuracy of M, this is unlikely to be met.
%   
%   See also CRRECTMAP, EVAL.

%   Copyright 1998 by Toby Driscoll.
%   $Id: evalinv.m 7 1998-05-10 04:37:19Z tad $

param = parameters(M.diskmap);

if nargin < 3
  qdata = param.qdata;
  qdatar = M.rectqdata;
else 
  qdatar = tol;
  qdata = [];
end

p = polygon(M.diskmap);
w = vertex(p);
beta = angle(p) - 1;
pr = M.rectpolygon;
wr = vertex(pr);
betar = angle(pr) - 1;
cr = param.crossratio;
aff = param.affine;
affr = M.rectaffine;
Q = param.qlgraph;

wp = NaN*zp;
idx = logical(ones(size(zp)));
wp(idx) = crrmap(zp(idx),wr,betar,w,beta,cr,affr,aff,Q,qdatar,qdata);
