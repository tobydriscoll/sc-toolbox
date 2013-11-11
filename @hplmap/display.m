function out = display(M)
%DISPLAY Display parameters of a Schwarz-Christoffel half-plane map.

%   Copyright 1998 by Toby Driscoll.
%   $Id: display.m 118 2001-05-03 21:18:27Z driscoll $

p = polygon(M);
w = vertex(p);
alpha = angle(p);
z = M.prevertex;
c = M.constant;

if length(z) < length(w)
  z = [z(:);Inf];
end

L = {' '; '  hplmap object:'; ' '};
L{4} = '      vertex               alpha         prevertex       ';
L{5} = ' --------------------------------------------------------';
u = real(w);
v = imag(w);
fmt = ' %8.5f %c %7.5fi     %8.5f    %20.12e';
for j = 1:length(w)
  if v(j) < 0
    s = '-';
  else
    s = '+';
  end
  L{end+1} = sprintf(fmt,u(j),s,abs(v(j)),alpha(j),z(j));
end
L{end+1} = ' ';
if imag(c) < 0
  s = '-';
else
  s = '+';
end
L{end+1} = sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c)));
L{end+1} = ' ';
L{end+1} = sprintf('  Apparent accuracy is %.2e',M.accuracy);
L{end+1} = ' ';


if nargout==0
  fprintf('%s\n',L{:})
else
  out = L;
end