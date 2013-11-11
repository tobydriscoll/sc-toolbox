function v = parameters(M)
%PARAMETERS Return a structure of the Schwarz-Christoffel map parameters.

%   Copyright 1998 by Toby Driscoll.
%   $Id: parameters.m 7 1998-05-10 04:37:19Z tad $

v.crossratio = M.crossratio;
v.affine = M.affine;
v.center = M.center{1};
v.qlgraph = M.qlgraph;
v.original = M.original;
v.qdata = M.qdata;

% We compute prevertices, even though they aren't used.
wcfix = M.center{2};
mt = wcfix(2:5);
z = embedding(M,wcfix(1));
v.prevertex = (mt(4)*z - mt(2)) ./ (-mt(3)*z + mt(1));
