function Md = diskmap(Ms)
%   Convert strip map to disk map.

%   Copyright 1998 by Toby Driscoll.
%   $Id: diskmap.m 43 1998-07-01 19:12:31Z tad $

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
zh = exp(pi*z);

% Map -Inf correctly
idx = find(isinf(z) & (z < 0));
zh(idx) = 0;

% Map Inf correctly
zh(n) = Inf;

% Transform prevertices to disk
A = moebius(zh(n-2:n),[-1 -i 1]);
zd = sign(A(zh));
zd(n) = 1;

% Create new map
Md = diskmap(polygon(w,beta+1),zd);
