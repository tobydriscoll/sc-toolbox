function p = subsasgn(p,S,data)
%   Allows assignment of individual vertices.

%   Copyright 1998 by Toby Driscoll.
%   $Id: subsasgn.m 7 1998-05-10 04:37:19Z tad $

name = fieldnames(p);

% Single index reference p(idx)
if length(S) == 1 & strcmp(S.type,'()') & length(S.subs) == 1
  p.vertex(S.subs{1}) = data;
else
  % Field reference
  if strcmp(S(1).type,'.')
    field = S(1).subs;
    idx = strmatch(field,name,'exact');
    if isempty(idx)
      error(sprintf('Unrecognized field name ''%s''.',field));
    end
    if length(S) > 1 & strcmp(S(2).type,'()') & length(S(2).subs) == 1
      idx = S(2).subs{1};
    else
      idx = ':';
    end
    eval(sprintf('p.%s(idx) = data;',field))
  end
end  