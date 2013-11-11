function corner = corners(M)
%CORNERS Indices of rectangle/generalized quadrilateral corners.

%   Copyright 1998 by Toby Driscoll.
%   $Id: corners.m 7 1998-05-10 04:37:19Z tad $

z = M.prevertex;
tol = 4*eps;

% Find extent of rectangle
K = max(real(z));
Kp = max(imag(z));

% First corner is K + 0i
dif = repmat(z,1,4) - repmat([K K+i*Kp -K+i*Kp -K],length(z),1);
[tmp,corner] = min(abs(dif));

corner = corner(:);
