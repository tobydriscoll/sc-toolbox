function [h,r,theta] = plot(f,varargin)
%PLOT Visualize a Schwarz-Christoffel Riemann surface map.
%   PLOT(F) plots the polygon associated with the Schwarz-Christoffel
%   Riemann surface map F and the images of ten evenly spaced circles
%   and radii under the S-C transformation.
%   
%   PLOT(F,NR,NTHETA) plots the images of NR circles and NTHETA radii.
%   
%   PLOT(F,R,THETA) plots the circles at radii given by the entries of R
%   and radii at the angles specified in THETA.
%   
%   PLOT(F,TOL) or PLOT(F,NR,NTHETA,TOL) or PLOT(F,R,THETA,TOL) computes
%   the map with accuracy roughly TOL. Normally TOL defaults to 1e-4 or
%   the accuracy of F, whichever is greater.
%   
%   See also RIESURFMAP, EVAL.

%   Copyright 2002 by Toby Driscoll.
%   $Id: plot.m 298 2009-09-15 14:36:37Z driscoll $

p = polygon(f);
w = vertex(p);
beta = angle(p) - 1;
z = f.prevertex;
zb = f.prebranch;
c = f.constant;

if nargin == 1
  [a1,a2,a3] = rsplot(w,beta,z,zb,c);
elseif length(varargin) == 1
  % Tolerance given only
  [a1,a2,a3] = rsplot(w,beta,z,zb,c,10,10,ceil(-log10(varargin{1})));
elseif length(varargin) == 2
  % R, theta given only
  [a1,a2,a3] = rsplot(w,beta,z,zb,c,varargin{1},varargin{2});
else
  % All given
  nqpts = ceil(-log10(varargin{3}));
  [a1,a2,a3] = rsplot(w,beta,z,zb,c,varargin{1},varargin{2},nqpts);
end

if nargout > 0
  h = a1;
  r = a2;
  theta = a3;
end

