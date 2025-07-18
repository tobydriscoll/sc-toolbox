function Mh = hplmap(Ms)
%   Convert strip map to half-plane map.

%   Copyright 1998 by Toby Driscoll.
%   $Id: hplmap.m 7 1998-05-10 04:37:19Z tad $

p = polygon(Ms);
w = vertex(p);
beta = angle(p) - 1;
z = Ms.prevertex;
n = length(z);

% Find index of vertex at Inf
idx = find(isinf(z) & (z > 0));

% Put that vertex last
renum = [idx+1:n 1:idx];
w = w(renum);
beta = beta(renum);
z = z(renum);

% Map prevertices to real axis
zh = real(exp(pi*z));

% Map -Inf correctly
idx = isinf(z) & (z < 0);
zh(idx) = 0;

% Put finite ones inside [-1,1]
zh = (zh-zh(1)) * 2/(zh(n-1)-zh(1)) - 1;

% Map Inf correctly
zh(n) = Inf;

% Create new map
Mh = hplmap(polygon(w,beta+1),zh);
