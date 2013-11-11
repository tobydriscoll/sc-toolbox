function zr = rectangle(M)
%RECTANGLE Return the corners of the rectangle in the fundamental domain.

%   Copyright 1998 by Toby Driscoll.
%   $Id: rectangle.m 7 1998-05-10 04:37:19Z tad $

zr = prevertex(M);
zr = zr(corners(M));
