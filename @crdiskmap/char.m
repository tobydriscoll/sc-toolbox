function out = char(f)
%CHAR   Pretty-print a Schwarz-Christoffel crossratio disk map.

%   Copyright 1998-2001 by Toby Driscoll.
%   $Id: char.m 162 2001-07-20 14:33:00Z driscoll $

p = polygon(f);
w = vertex(p);
beta = angle(p)-1;
cr = f.crossratio;
Q = f.qlgraph;

L = cell(2,1);
L{1}='   Quadrilateral vertices        Prevertex crossratio  ';
L{2}=' ------------------------------------------------------';
for j = 1:length(cr)
  L{end+1}=sprintf('      %2i  %2i  %2i  %2i                   %8.5f',...
      Q.qlvert(:,j),cr(j));
end
wc = center(f);
if imag(wc) < 0
  s = '-';
else
  s = '+';
end
L{end+1} = ' ';
L{end+1} = sprintf('   Conformal center at %.4f %c %.4fi',real(wc),s,abs(imag(wc)));
L{end+1} = sprintf('   Apparent accuracy is %.2e',f.accuracy);
L{end+1} = ' ';

out = L;
