function d = diam(p)
%DIAM    Diameter of a polygon.
%
%   DIAM(P) returns max_{j,k} |P(j)-P(k)|. This may be infinite.

%   Copyright 2002 by Toby Driscoll.
%   $Id: diam.m 194 2002-09-10 19:10:41Z driscoll $

w = vertex(p);
[w1,w2] = meshgrid(w);
d = max( max( abs(w1-w2) ) );
