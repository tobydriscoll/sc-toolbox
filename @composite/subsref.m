function wp = subsref(f,S)
%SUBSREF Evaluate composite map by subscript notation.
%   F(ZP) is a synonym for EVAL(F,ZP).
%   
%   See also EVAL, COMPOSITE.

%   Copyright 2001 by Toby Driscoll.
%   $Id: subsref.m 154 2001-07-20 13:52:46Z driscoll $

if length(S) == 1 & strcmp(S.type,'()')
  wp = eval(f,S.subs{1});
elseif length(S) == 1 & strcmp(S.type,'{}')
  idx = S.subs{1};
  if length(idx)==1 & ~ischar(idx)
    wp = f.maps{idx};
  else
    wp = f.maps(idx);
  end
end

  