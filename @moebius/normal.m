function M = normal(M)
%NORMAL Normalize a Moebius transformation.
%   NORMAL(M), where M is a moebius object, divides the coefficients of M
%   by the constant term in the demoninator, or the coefficient of the
%   linear term in the denominator if the constant is zero.
%
%   Copyright (c) 2004 by Toby Driscoll.
%   $Id: normal.m 273 2004-01-23 17:25:17Z driscoll $

c = M.coeff(3);
if c==0, c = M.coeff(4); end
M.coeff = M.coeff/c;

