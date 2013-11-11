function M = uminus(M)
%UMINUS Negate a Moebius transformation.
%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: uminus.m 44 1998-07-01 20:21:54Z tad $

M.coeff(1:2) = -M.coeff(1:2);
