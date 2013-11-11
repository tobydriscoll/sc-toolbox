function M = mrdivide(M1,M2)
%MRDIVIDE Divide Moebius map by a scalar, or reciprocate it.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: mrdivide.m 44 1998-07-01 20:21:54Z tad $

if isa(M1,'double') & length(M1)==1 
  % Exchange numerator and denominator
  M2.coeff = M2.coeff([3 4 1 2]);
  % Multiply by scalar
  M = M1*M2;
elseif isa(M2,'double') & length(M2)==1
  C = M1.coeff;
  C(3:4) = M2*C(3:4);
  M = moebius;
  M.coeff = C;
else
  error('Division not defined for these operands.')
end
