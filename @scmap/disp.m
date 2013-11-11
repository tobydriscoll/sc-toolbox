function disp(f)

% $Id: disp.m 156 2001-07-20 14:03:14Z driscoll $ 

s = char(f);
if isstr(s)
  disp(s)
elseif iscell(s)
  fprintf('\n  SC %s:\n\n',class(f));
  for n = 1:length(s)
    disp(s{n})
  end
end
