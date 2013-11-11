function Mi = scmapinv(M)
%SCMAPINV Inverse of a Schwarz-Christoffel map.
%   SCMAPINV(M) returns a dummy object that represents the inverse of
%   the SC map M. The only things possible with this object are to EVAL
%   it or INV it back to the SC map.
%   
%   See also SCMAPINV/EVAL, SCMAPINV/SUBSREF, SCMAPINV/INV.

%   Copyright 1998 by Toby Driscoll.
%   $Id: scmapinv.m 155 2001-07-20 13:53:44Z driscoll $

if nargin==0
  Mi.themap = [];
else
  Mi.themap = M;
end

Mi = class(Mi,'scmapinv');
