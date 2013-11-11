function out = display(M)
%DISPLAY Display parameters of a Schwarz-Christoffel CR rectified map.

%   Copyright 1998 by Toby Driscoll.
%   $Id: display.m 118 2001-05-03 21:18:27Z driscoll $

p = polygon(M);
w = vertex(p);
wr = vertex(M.rectpolygon);
alphar = angle(M.rectpolygon);

L = {' '; '  crrectmap object:'; ' '};
L{4} = '      vertex          rectified angle        prevertex       ';
L{5} = ' ------------------------------------------------------------';

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
L{end+1} = sprintf('  Apparent accuracy is %.2e',accuracy(M));
L{end+1} = ' ';


if nargout==0
  fprintf('%s\n',L{:})
else
  out = L;
end
