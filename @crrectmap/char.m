function out = char(f)
%CHAR   Display parameters of a Schwarz-Christoffel CR rectified map.

%   Copyright 1998-2001 by Toby Driscoll.
%   $Id: char.m 162 2001-07-20 14:33:00Z driscoll $

p = polygon(f);
w = vertex(p);
wr = vertex(f.rectpolygon);
alphar = angle(f.rectpolygon);

L = cell(2,1);
L{1} = '      vertex          rectified angle        prevertex       ';
L{2} = ' ------------------------------------------------------------';

u = real(w);
v = imag(w);
ur = real(wr);
vr = imag(wr);
sgn = ['-','+','+'];
s = sgn(sign(v)+2);
sr = sgn(sign(vr)+2);
fmt = ' %8.5f %c %7.5fi       %3.1f pi       %8.5f %c %7.5fi';
for j = 1:length(w)
  L{end+1} = sprintf(fmt,...
      u(j),s(j),abs(v(j)),alphar(j),ur(j),sr(j),abs(vr(j)));
end
L{end+1} = ' ';
L{end+1} = sprintf('  Apparent accuracy is %.2e',accuracy(f));
L{end+1} = ' ';

out = L;
