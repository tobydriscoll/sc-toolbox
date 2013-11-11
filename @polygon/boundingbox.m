function box = boundingbox(p)
%BOUNDINGBOX Smallest box that contains the polygon.
%   BOUNDINGBOX(P) returns the smallest box (in AXIS format) that contains
%   the polygon P. If P is unbounded, all the entries will be infinite.
%
%   See also POLYGON/DIAM.
%
%   Copyright 2003 by Toby Driscoll.
%   $Id: boundingbox.m 266 2003-04-25 18:46:31Z driscoll $

if ~isinf(p)
  z = vertex(p);
  box = [ min(real(z)) max(real(z)) min(imag(z)) max(imag(z)) ];
else
  % We might find some finite bounds. But is there any application for this?
  box = [ -Inf Inf -Inf Inf ];
end