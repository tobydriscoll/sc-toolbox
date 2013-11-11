function out = display(M)
%DISPLAY Display parameters of a Schwarz-Christoffel exterior map.

%   Copyright 1998 by Toby Driscoll.
%   $Id: display.m 118 2001-05-03 21:18:27Z driscoll $

p = polygon(M);
w = vertex(p);
alpha = angle(p);
n = length(w);
z = flipud(M.prevertex);
c = M.constant;

L = { ' '; '  extermap object:'; ' ' };
L{4}='      vertex              alpha         prevertex             arg/pi';
L{5}=' -----------------------------------------------------------------------';
u = real(w);
v = imag(w);
x = real(z);
y = imag(z);
ang = angle(z)/pi;
ang(ang<=0) = ang(ang<=0) + 2;
fmt = ' %8.5f %c %7.5fi    %8.5f   %8.5f %c %7.5fi    %14.12f';
for j = 1:length(w)
  if v(j) < 0
    s1 = '-';
  else
    s1 = '+';
  end
  if y(j) < 0
    s2 = '-';
  else
    s2 = '+';
  end
  L{end+1}=sprintf(fmt,u(j),s1,abs(v(j)),alpha(j),x(j),s2,abs(y(j)),ang(j));
end

L{end+1} = ' ';
if imag(c) < 0
  s = '-';
else
  s = '+';
end
L{end+1} = sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c)));
L{end+1} = sprintf('  Logarithmic capacity = %.8g',abs(c));
L{end+1} = ' ';
L{end+1} = sprintf('  Apparent accuracy is %.2e',M.accuracy);
L{end+1} = ' ';


if nargout==0
  fprintf('%s\n',L{:})
else
  out = L;
end
