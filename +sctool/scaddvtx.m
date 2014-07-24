function [wn,betan] = scaddvtx(w,beta,pos,window)
%SCADDVTX Add a vertex to a polygon.
%   [WN,BETAN] = SCADDVTX(W,BETA,POS) adds a new vertex to the polygon
%   described by W and BETA immediately after vertex POS.  If
%   W(POS:POS+1) are finite, the new vertex is at the midpoint of an
%   edge; otherwise, the new vertex is a reasonable distance from its
%   finite neighbor.
%       
%   See also SCFIX.
    
%   Copyright 1998 by Toby Driscoll.
%   $Id: scaddvtx.m 298 2009-09-15 14:36:37Z driscoll $

if nargin < 4, window = [-Inf Inf -Inf Inf]; end
w = w(:);
beta = beta(:);
n = length(w);
if ~pos, pos=n; end
pos1 = rem(pos,n)+1;
if ~any(isinf(w([pos,pos1])))	% easy case
  new = mean(w([pos,pos1]));
else					% messy case
  % Find a pair of adjacent finite vertices as a basis.
  base = min(find(~isinf(w) & ~isinf(w([2:n,1]))));
  ang(base) = angle(w(rem(base,n)+1)-w(base));
  
  % Determine absolute angle of side pos->pos1.
  for j = [base+1:n,1:base-1]
    ang(j) = ang(rem(j-2+n,n)+1)-pi*beta(j);
    if j==pos, break, end
  end
  
  % Find a nice side length.
  len = abs(w([2:n,1])-w);
  avglen = mean(len(~isinf(len)));

  if isinf(w(pos))
    base = w(pos1);
    dir = exp(i*(ang(pos)+pi));
  else
    base = w(pos);
    dir = exp(i*(ang(pos)));
  end
  
  % Restrict point to a window (to help out graphics).
  new = base + avglen*dir;  
  while real(new) < window(1) | real(new) > window(2) | ...
	imag(new) < window(3) | imag(new) > window(4)
    avglen = avglen / 2;
    new = base + avglen*dir;
  end
    
end

wn = [w(1:pos);new;w(pos+1:n)];
betan = [beta(1:pos);0;beta(pos+1:n)];
