function [h,r,theta] = plot(M,varargin)
%PLOT Visualize a Schwarz-Christoffel disk map.
%   PLOT(M) plots the polygon associated with the Schwarz-Christoffel
%   disk map M and the images of ten evenly spaced circles and radii
%   under the S-C transformation.
%   
%   PLOT(M,NR,NTHETA) plots the images of NR circles and NTHETA radii.
%   
%   PLOT(M,R,THETA) plots the circles at radii given by the entries of R
%   and radii at the angles specified in THETA.
%   
%   PLOT(M,TOL) or PLOT(M,NR,NTHETA,TOL) or PLOT(M,R,THETA,TOL)
%   computes the map with accuracy roughly TOL. Normally TOL defaults to
%   1e-4 or the accuracy of M, whichever is greater.
%   
%   See also DISKMAP, EVAL.

%   Copyright 1998 by Toby Driscoll.
%   $Id: plot.m 7 1998-05-10 04:37:19Z tad $

p = polygon(M);
w = vertex(p);
beta = angle(p) - 1;
z = M.prevertex;
c = M.constant;

if nargin == 1
  [a1,a2,a3] = dplot(w,beta,z,c);
elseif length(varargin) == 1
  % Tolerance given only
  [a1,a2,a3] = dplot(w,beta,z,c,10,10,ceil(-log10(varargin{1})));
elseif length(varargin) == 2
  % R, theta given only
  [a1,a2,a3] = dplot(w,beta,z,c,varargin{1},varargin{2});
else
  % All given
  nqpts = ceil(-log10(varargin{3}));
  [a1,a2,a3] = dplot(w,beta,z,c,varargin{1},varargin{2},nqpts);
end

if nargout > 0
  h = a1;
  r = a2;
  theta = a3;
end

