function out = display(M)
%DISPLAY Display parameters of a Schwarz-Christoffel rectangle map.

%   Copyright 1998 by Toby Driscoll.
%   $Id: display.m 118 2001-05-03 21:18:27Z driscoll $

p = polygon(M);
w = vertex(p);
alpha = angle(p);
z = M.prevertex;
c = M.constant;

n = length(w);
% Deduce corner locations
left = abs(real(z)-min(real(z))) < eps;
right = abs(real(z)-max(real(z))) < eps;
top = abs(imag(z)-max(imag(z))) < eps;
bot = abs(imag(z)-min(imag(z))) < eps;
corners = find(left+right+top+bot - 1);
c1 = find(abs(z-max(real(z))) < eps);
offset = find(corners==c1);
corners = corners([offset:4,1:offset-1]);
rect = z(corners);

L = {' '; '  rectmap object:'; ' '};
L{4}=' cnr      vertex              alpha               prevertex       ';
L{5}=' ------------------------------------------------------------------------';

u = real(w);
v = imag(w);
for j = 1:length(w)
  if v(j) < 0
    s = '-';
  else
    s = '+';
  end
  cnr = find(j==corners);
  if isempty(cnr)
    cstr = '    ';
  else
    cstr = sprintf('  %i ',cnr);
  end
  if ~imag(z(j))
    L{end+1}=sprintf('%s %8.5f %c %7.5fi    %8.5f   %16.8e',...
	cstr,u(j),s,abs(v(j)),alpha(j),z(j));
  else
    L{end+1}=sprintf('%s %8.5f %c %7.5fi    %8.5f   %16.8e + %14.8ei',...
      cstr,u(j),s,abs(v(j)),alpha(j),real(z(j)),imag(z(j)));
  end    
end

L{end+1} = ' ';
if imag(c) < 0
  s = '-';
else
  s = '+';
end
L{end+1} = sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c)));
L{end+1} = sprintf('  Conformal modulus = %.8g',imag(rect(2))/rect(1)/2);
L{end+1} = ' ';
L{end+1} = sprintf('  Apparent accuracy is %.2e',M.accuracy);
L{end+1} = ' ';


if nargout==0
  fprintf('%s\n',L{:})
else
  out = L;
end