function q = truncate(p)
%TRUNCATE Truncate an unbounded polygon.
%   Q = TRUNCATE(P) returns a polygon whose finite vertices are the same
%   as those of P and whose infinite vertices have been replaced by
%   several finite ones. The new vertices are chosen by using a
%   square "cookie cutter" on P.

%   Copyright 2002 by Toby Driscoll. All rights reserved.
%   $Id: truncate.m 229 2003-01-09 14:45:20Z driscoll $

w = vertex(p);
n = length(w);

if ~any(isinf(w)), q=p; return, end

% Put the finite vertices in a box.
wf = w(~isinf(w));
xbound = [ min(real(wf)); max(real(wf)) ];
ybound = [ min(imag(wf)); max(imag(wf)) ];
delta = max( diff(xbound), diff(ybound) );
xrect = mean(xbound) + [-1;1]*delta;
yrect = mean(ybound) + [-1;1]*delta;
zrect = xrect([1 1 2 2]') + 1i*yrect([2 1 1 2]');

% Find intersections of the box with the unbounded sides.
[hit,loc] = intersect(p,[zrect zrect([2:4 1])]);

% Carry over the finite vertices, inserting to substitute for the
% infinite ones.
atinf = find(isinf(w));
v = w(1:atinf(1)-1);
for k = 1:length(atinf)
  M = atinf(k); M1 = mod(M-2,n) + 1;
  % Find where the adjacent sides hit the rectangle.
  sp = find( hit(:,M1) );  rp = loc(sp,M1);
  sn = find( hit(:,M) );   rn = loc(sn,M);
  % Include the rectangle corners that are "in between".
  dt = mod( angle(rn/rp), 2*pi );
  dr = mod( angle(zrect/rp), 2*pi );
  [dr,idx] = sort(dr);
  use = dr<dt;
  v = [ v; rp; zrect(idx(use)); rn ];
  if k < length(atinf)
    v = [ v; w(M+1:atinf(k+1)-1) ];
  else
    v = [ v; w(M+1:end) ];
  end
end

q = polygon(v);