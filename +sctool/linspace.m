function y = linspace(d1, d2, n)
%LINSPACE Linearly spaced vector.
%   LINSPACE(x1, x2) generates a row vector of 100 linearly
%   equally spaced points between x1 and x2.
%
%   LINSPACE(x1, x2, N) generates N points between x1 and x2.
%
%   See also LOGSPACE, :.
%
%   This version is modified from that provided in MATLAB 5.2,
%   because that version fails when the first two arguments are
%   complex. 
%   
%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: linspace.m 2 1998-05-10 02:57:52Z tad $

if nargin == 2
    n = 100;
end
if n~=1
  % Original: y = d1:(d2-d1)/(n-1):d2;
  y = d1 + (d2-d1)*(0:1/(n-1):1);
else
  y = d2;
end
