function M = minus(M1,M2)
%MINUS  Subtract a scalar from a Moebius map.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: minus.m 44 1998-07-01 20:21:54Z tad $

M = M1 + (-M2);
