function rgn = dscpolygons(p1,p0)
%DSCPOLYGONS Construct doubly connected regions.
%   DSCPOLYGONS(P1,P0) constructs a region which is bounded by an inner
%   polygon P1 and an outer polygon P0.
%
%   See also POLYGON.
%   Written by Alfa Heryudono, 2003.
  
superiorto('double');

if nargin ~= 2
    error('Wrong input')
end

%   Check if p1 is "inside" p0. In other words, check for intesections between sides of the
%   polygon p1 and the polygon p0.
%   Not implemented yet.

%   Check if p1 and p0 are empty
%   Not implemented yet

%   Create a region bounded by p1 and p0
rgn.p1 = p1;
rgn.p0 = p0;
rgn = class(rgn,'dscpolygons');


