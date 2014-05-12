function varargout = get(rgn)
%GET Get outer and inner polygons of a doubly connected region.

%   [P1,P0] = GET(RGN) returns the values of the inner polygon P1 and the
%   values of the outer polygon P0

%   Written by Alfa Heryudono, 2003.

varargout{1} = rgn.p1;
varargout{2} = rgn.p0;