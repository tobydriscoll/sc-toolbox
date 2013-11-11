function d = crpsdist(segment,pts)
%CRPSDIST Distance from point(s) to a line segment.
%   CRPSDIST(SEG,PTS) returns a vector the size of PTS indicating the
%   distance from each entry to the line segment described by SEG.

%   Copyright 1997 by Toby Driscoll. Last updated 04/29/97.

if isempty(pts)
  d = [];
  return
end

d = Inf*pts;

% Rotate to make segment equal [0,xmax].
pts = (pts-segment(1))/sign(diff(segment));
xmax = abs(diff(segment));

% Some points are closest to segment's interior.
mask = (real(pts)>=0) & (real(pts)<=xmax);
d(mask) = abs(imag(pts(mask)));

% The others are closest to an endpoint.
mask = real(pts) < 0;
d(mask) = abs(pts(mask));
mask = real(pts) > xmax;
d(mask) = abs(pts(mask)-xmax);
