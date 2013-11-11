function [x,y] = vertex(p)
%VERTEX Vertices of a polygon.
%   VERTEX(P) returns the vertices of polygon P as a complex vector.
%   
%   [X,Y] = VERTEX(P) returns the vertices as two real vectors.
%   
%   See also POLYGON.

%   Copyright 1998 by Toby Driscoll.
%   $Id: vertex.m 7 1998-05-10 04:37:19Z tad $

x = p.vertex;
if nargout == 2
  y = imag(x);
  x = real(x);
end
