function Md = diskmap(M)
%DISKMAP Convert to diskmap object.
%   DISKMAP(M), where M is a crossratio diskmap object, returns the
%   equivalent map as a diskmap object. This is done by computing the
%   prevertices which yield the correct conformal center and have the
%   last prevertex equal to 1. If the crossratios of M exceed about 20,
%   accuracy will probably be lost in some parts of the polygon.

%   Copyright 1998 by Toby Driscoll.
%   $Id: diskmap.m 7 1998-05-10 04:37:19Z tad $

cr = M.crossratio;
wcfix = M.center{2};

% Prevertices in the embedding described in wcfix
z = crembed(cr,M.qlgraph,wcfix(1));

% Transform to make conformal center correct
mt = wcfix(2:5);
z = (-mt(4)*z + mt(2))./(mt(3)*z - mt(1));

Md = diskmap(polygon(M),z);
