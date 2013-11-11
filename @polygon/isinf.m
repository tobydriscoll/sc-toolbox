function s = isinf(p)
%ISINF   Determine whether a polygon is unbounded.
%
%   ISINF(P) returns logical 1 when P is unbounded and 0 otherwise.

%   Copyright 2002 by Toby Driscoll.
%   $Id: isinf.m 190 2002-09-10 17:23:26Z driscoll $

s = any(isinf(vertex(p)));