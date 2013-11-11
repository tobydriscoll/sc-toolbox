function handles = plot(p,varargin)
%PLOT Plot a polygon.
%   PLOT(P) plots the polygon P (draws the sides).
%   
%   PLOT(P,'num') or PLOT(P,'lab') also plots dots for the vertices and
%   numeric labels. For infinite vertices, two numeric labels are
%   printed.
%   
%   H = PLOT(P) returns a vector of handles to the sides drawn. H =
%   PLOT(P,'num') returns a structure; H.side is a vector of side
%   handles, and H.label is a 2-column array of vertex dot and label
%   handles (each vertex has two graphical objects).

%   Copyright 1998 by Toby Driscoll.
%   $Id: plot.m 138 2001-05-15 14:16:46Z driscoll $

w = vertex(p);
beta = angle(p) - 1;

if isempty(w)
  h = [];
elseif nargin > 1
  [eh,lh] = plotpoly(w,beta,1);
  h.side = eh;
  h.label = lh;
else
  h = plotpoly(w,beta);
end

if nargout > 0
  handles = h;
end

