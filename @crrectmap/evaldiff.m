function fp = evaldiff(M,zp,tol)
%EVALDIFF Derivative of Schwarz-Christoffel crossratio disk map at points.
%   EVALDIFF(M,ZP) computes the derivative of the Schwarz-Christoffel
%   disk map M at the points ZP.
%   
%   See also CRDISKMAP, EVAL.

%   Copyright 1998 by Toby Driscoll.
%   $Id: evaldiff.m 87 1999-11-19 21:51:52Z tad $

% Special thanks to Stefano Costa, who actually made this work.


param = parameters(M.diskmap);

if nargin < 3
  qdata = param.qdata;
  qdatar = M.rectqdata;
else 
  qdata = tol;
  qdatar = [];
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

fp = crrderiv(zp,w,beta,wr,betar,cr,aff,affr,Q,qdata,qdatar);
