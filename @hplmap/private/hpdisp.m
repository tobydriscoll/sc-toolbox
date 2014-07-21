function hpdisp(w,beta,z,c)
%HPDISP Display results of Schwarz-Christoffel half-plane parameter problem.
%   HPDISP(W,BETA,Z,C) displays the results of HPPARAM in a pleasant
%   way.
%
%   See also HPPARAM, HPPLOT.

%   Copyright 1998 by Toby Driscoll.
%   $Id: hpdisp.m 298 2009-09-15 14:36:37Z driscoll $

if length(z) < length(w)
  z = [z(:);Inf];
end
disp(' ')
disp('      vertex [w]           beta          prevertex [z]   ')
disp(' --------------------------------------------------------')
u = real(w);
v = imag(w);
for j = 1:length(w)
  if v(j) < 0
    s = '-';
  else
    s = '+';
  end
  disp(sprintf(' %8.5f %c %7.5fi     %8.5f    %20.12e',...
      u(j),s,abs(v(j)),beta(j),z(j)));
end
disp(' ')
if imag(c) < 0
  s = '-';
else
  s = '+';
end
disp(sprintf('  c = %.8g %c %.8gi\n',real(c),s,abs(imag(c))))

