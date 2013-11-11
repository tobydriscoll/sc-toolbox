function [h,re,im] = plot(M,varargin)
%PLOT Visualize a Schwarz-Christoffel crossratio disk map.
%   PLOT(M) plots the polygon associated with the Schwarz-Christoffel
%   crossratio rectified map M and the images of ten evenly spaced
%   vertical and horizontal line segments the S-C transformation.
%   
%   PLOT(M,NRE,NIM) plots the images of NRE vertical and NIM horizontal
%   line segments.
%   
%   PLOT(M,RE,IM) plots the vertical line segments at abscissae given by
%   the entries of RE and horizontal line segments at the ordinates
%   specified in IM.
%   
%   PLOT(M,TOL) or PLOT(M,NRE,NIM,TOL) or PLOT(M,RE,IM,TOL) computes the
%   map with accuracy roughly TOL. Normally TOL defaults to 1e-4 or the
%   accuracy of M, whichever is greater.
%   
%   See also CRRECTMAP, EVAL.

%   Copyright 1998 by Toby Driscoll.
%   $Id: plot.m 7 1998-05-10 04:37:19Z tad $

param = parameters(M.diskmap);
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

if nargin == 1
  [a1,a2,a3] = crrplot(w,beta,wr,betar,cr,aff,affr,Q);
elseif length(varargin) == 1
  % Tolerance given only
  [a1,a2,a3] = crrplot(w,beta,wr,betar,cr,aff,affr,Q,...
      10,10,ceil(-log10(varargin{1})));
elseif length(varargin) == 2
  % re,im given only
  [a1,a2,a3] = crrplot(w,beta,wr,betar,cr,aff,affr,Q,varargin{:});
else
  % All given
  nqpts = ceil(-log10(varargin{3}));
  [a1,a2,a3] = crrplot(w,beta,wr,betar,cr,aff,affr,Q,...
      varargin{1},varargin{2},nqpts);
end

if nargout > 0
  h = a1;
  re = a2;
  im = a3;
end
