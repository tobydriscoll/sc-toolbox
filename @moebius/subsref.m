function w = subsref(M,S)
%SUBSREF Evaluate or compose maps.
%   M(Z), where M is a Moebius transformation and Z is a vector of
%   points, returns the image of the points in Z. This just a synonym
%   for EVAL(M,Z). 
%   
%   M1(M2), where M1 and M2 are both Moebius transformations, returns
%   a new Moebius transformation that is their composition.
%   
%   See also EVAL, MOEBIUS.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: subsref.m 44 1998-07-01 20:21:54Z tad $

if length(S) == 1 & strcmp(S.type,'()')
  if isa(S.subs{1},'double')
    w = eval(M,S.subs{1});
  elseif isa(S.subs{1},'moebius')
    C = M.coeff;
    D = S.subs{1}.coeff;
    w = moebius;
    w.coeff = [C(1:2)*D([3 1]).' C(1:2)*D([4 2]).' ...
          C(3:4)*D([3 1]).' C(3:4)*D([4 2]).'];
  else
    error('Syntax not supported.')
  end
else
  error('Only syntax for MOEBIUS is a single parenthesized subscript.')
end
