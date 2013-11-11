function v = parameters(M)
%PARAMETERS Return a structure of the Schwarz-Christoffel map parameters.

%   Copyright 1998 by Toby Driscoll.
%   $Id: parameters.m 7 1998-05-10 04:37:19Z tad $

v.diskmap = M.diskmap;
v.rectpolygon = M.rectpolygon;
v.rectaffine = M.rectaffine;
v.prevertex = vertex(M.rectpolygon);
