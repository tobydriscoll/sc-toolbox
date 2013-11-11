function out = char(f)
%CHAR   Pretty-print a Schwarz-Christoffel strip map.

%   Copyright 1998-2001 by Toby Driscoll.
%   $Id: char.m 162 2001-07-20 14:33:00Z driscoll $

p = polygon(f);
w = vertex(p);
alpha = angle(p);
z = f.prevertex;
c = f.constant;

L = cell(2,1);
L{1} = '      vertex               alpha           prevertex         ';
L{2}=' ------------------------------------------------------------';
u = real(w);
v = imag(w);
for j = 1:length(w)
  if v(j) < 0
    s = '-';
  else
    s = '+';
  end
  if ~imag(z(j))
    L{end+1}=sprintf(' %8.5f %c %7.5fi     %8.5f    %20.12e',...
	u(j),s,abs(v(j)),alpha(j),z(j));
  else
    L{end+1}=sprintf(' %8.5f %c %7.5fi     %8.5f    %20.12e + i',...
      u(j),s,abs(v(j)),alpha(j),z(j));
  end    
end

L{end+1} = ' ';
if imag(c) < 0
  s = '-';
else
  s = '+';
end
L{end+1} = sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c)));
L{end+1} = sprintf('  Apparent accuracy is %.2e',f.accuracy);
L{end+1} = ' ';

out = L;
