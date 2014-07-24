function rdisp(w,beta,z,c,L)
%RDISP Display results of Schwarz-Christoffel rectangle parameter problem.
%   RDISP(W,BETA,RECT,Z,C) displays the results of RPARAM in a pleasant
%   way.
%
%   See also RPARAM, RPLOT.

%   Copyright 1998 by Toby Driscoll.
%   $Id: rdisp.m 298 2009-09-15 14:36:37Z driscoll $

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

disp(' ')
disp(' cnr      vertex [w]          beta                prevertex [z]   ')
disp(' ------------------------------------------------------------------------')
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
    disp(sprintf('%s %8.5f %c %7.5fi    %8.5f   %16.8e',...
      cstr,u(j),s,abs(v(j)),beta(j),z(j)));
  else
    disp(sprintf('%s %8.5f %c %7.5fi    %8.5f   %16.8e + %14.8ei',...
      cstr,u(j),s,abs(v(j)),beta(j),real(z(j)),imag(z(j))));
  end    
end
disp(' ')
if imag(c) < 0
  s = '-';
else
  s = '+';
end
disp(sprintf('  c = %.8g %c %.8gi',real(c),s,abs(imag(c))))
disp(sprintf('\n  Conformal modulus = %.8g\n',imag(rect(2))/rect(1)/2));

