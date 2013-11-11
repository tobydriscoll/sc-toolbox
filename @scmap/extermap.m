function M = extermap(M)
%EXTERMAP Convert generic Schwarz-Christoffel map object to exterior map.
%   EXTERMAP(M) creates a extermap object based on the polygon and
%   options contained in M.
%   
%   See the EXTERMAP class documentation.

%   Copyright 1998 by Toby Driscoll.
%   $Id: extermap.m 7 1998-05-10 04:37:19Z tad $

M = extermap(M.polygon,M.options);
