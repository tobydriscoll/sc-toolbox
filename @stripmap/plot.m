function [h,re,im] = plot(M,varargin)
%PLOT Visualize a Schwarz-Christoffel strip map.
%   PLOT(M) plots the polygon associated with the Schwarz-Christoffel
%   strip map M and the images of ten evenly spaced vertical line
%   segments and horizontal lines under the S-C transformation.
%   
%   PLOT(M,NRE,NIM) plots the images of NRE vertical line segments and
%   NIM horizontal lines.
%   
%   PLOT(M,RE,IM) plots the vertical line segments at abscissae given by
%   the entries of RE and horizontal lines at the ordinates specified in
%   IM.
%   
%   PLOT(M,TOL) or PLOT(M,NRE,NIM,TOL) or PLOT(M,RE,IM,TOL) computes the
%   map with accuracy roughly TOL. Normally TOL defaults to 1e-4 or the
%   accuracy of M, whichever is greater.
%   
%   See also STRIPMAP, EVAL.

%   Copyright 1998 by Toby Driscoll.
%   $Id: plot.m 7 1998-05-10 04:37:19Z tad $

p = polygon(M);
w = vertex(p);
beta = angle(p) - 1;
z = M.prevertex;
c = M.constant;

if nargin == 1
  [a1,a2,a3] = stplot(w,beta,z,c);
elseif length(varargin) == 1
  % Tolerance given only
  [a1,a2,a3] = stplot(w,beta,z,c,10,10,ceil(-log10(varargin{1})));
elseif length(varargin) == 2
  % RE,IM given only
  [a1,a2,a3] = stplot(w,beta,z,c,varargin{1},varargin{2});
else
  % All given
  nqpts = ceil(-log10(varargin{3}));
  [a1,a2,a3] = stplot(w,beta,z,c,varargin{1},varargin{2},nqpts);
end

if nargout > 0
  h = a1;
  re = a2;
  im = a3;
end

