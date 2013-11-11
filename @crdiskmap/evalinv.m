function zp = evalinv(M,wp,tol)
%EVALINV Invert Schwarz-Christoffel crossratio disk map at points.
%   EVALINV(M,WP) evaluates the inverse of the Schwarz-Christoffel map
%   M at the points WP in the polygon. The default tolerance of M is
%   used.
%   
%   EVALINV(M,WP,TOL) attempts to give an answer accurate to TOL. If TOL
%   is smaller than the accuracy of M, this is unlikely to be met.
%   
%   See also CRDISKMAP, CRDISKMAP/EVAL.

%   Copyright 1998 by Toby Driscoll.
%   $Id: evalinv.m 7 1998-05-10 04:37:19Z tad $

if nargin < 3
  % Default means use value stored in map object
  qdata = M.qdata;
  tol = M.accuracy;
else
  % An argument was supplied. Is it qdata or a tolerance?
  qdata = tol;
  if length(tol) > 1
    tol = 10^(-size(qdata,1));
  end
end
    
p = polygon(M);
w = vertex(p);
beta = angle(p) - 1;
cr = M.crossratio;
aff = M.affine;
wcfix = M.center{2};
Q = M.qlgraph;

zp = NaN*wp;
idx = logical(isinpoly(wp,p));

zp(idx) = crinvmap(wp(idx),w,beta,cr,aff,wcfix,Q,qdata,[0 tol]);
