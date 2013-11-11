function Md = diff(M)
%DIFF   Differentiated SC map object.
%   DIFF(M) returns an object formally representing the derivative of
%   the map M.

%   Copyright 1998 by Toby Driscoll.
%   $Id: diff.m 7 1998-05-10 04:37:19Z tad $

Md = scmapdiff(M);
