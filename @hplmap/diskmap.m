function M1 = diskmap(M)
%DISKMAP Convert Schwarz-Christoffel half-plane map to a map from the disk.

%   Copyright 1998 by Toby Driscoll.
%   $Id: diskmap.m 7 1998-05-10 04:37:19Z tad $

p = polygon(M);
[z1,c1] = hp2disk(vertex(p),angle(p)-1,M.prevertex,M.constant);
M1 = diskmap(p,scmapopt(M),z1,c1);


