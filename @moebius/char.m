function out = char(map)
%CHAR   Pretty-print a Moebius transformation.

%   Copyright (c) 1998-2001 by Toby Driscoll.
%   $Id: char.m 161 2001-07-20 14:32:59Z driscoll $

% Numerator
num = '';
a = map.coeff(1);
if a~=0
  if isreal(a)
    num = [num num2str(a,4)];
  else
    num = [num '(' num2str(a,4) ')'];
  end
end
a = map.coeff(2);
if a~=0
  if ~isempty(num)
    if real(a) >= 0
      num = [num ' + '];
    else
      num = [num ' - '];
      a = -a;
    end
  end
  if isreal(a) & a~=1
    num = [num num2str(a,4) '*'];
  else
    num = [num '(' num2str(a,4) ')*'];
  end
  num = [num 'z'];
end

% Denominator
den = '';
a = map.coeff(3);
if a~=0
  if isreal(a)
    den = [den num2str(a,4)];
  else
    den = [den '(' num2str(a,4) ')'];
  end
end
a = map.coeff(4);
if a~=0
  if ~isempty(den)
    if real(a) >= 0
      den = [den ' + '];
    else
      den = [den ' - '];
      a = -a;
    end
  end
  if isreal(a) & a~=1
    den = [den num2str(a,4) '*'];
  else
    den = [den '(' num2str(a,4) ')*'];
  end
  den = [den 'z'];
end

L = [length(num),length(den)];
D = (max(L)-L)/2;
num = [blanks(floor(D(1))) num blanks(ceil(D(1)))];
den = [blanks(floor(D(2))) den blanks(ceil(D(2)))];
fline = repmat('-',1,max(L));

out = sprintf('\n  %s\n  %s\n  %s\n\n',num,fline,den);

	  