function data = subsref(p,S)
%   Extract one or more vertices by index.

%   Copyright 1998 by Toby Driscoll.
%   $Id: subsref.m 7 1998-05-10 04:37:19Z tad $

name = fieldnames(p);

% Single index reference p(idx)
if length(S) == 1 & strcmp(S.type,'()') & length(S.subs) == 1
  data = p.vertex(S.subs{1});
  return
end

% Field reference
data = [];
if strcmp(S(1).type,'.')
  field = S(1).subs;
  idx = strmatch(field,name,'exact');
  if isempty(idx)
    error(sprintf('Unrecognized field name ''%s''.',field));
  end
  eval(sprintf('data = p.%s;',field'))
end

% Index on the reference
if length(S) > 1 & strcmp(S(2).type,'()') & length(S(2).subs) == 1
  data = data(S(2).subs{1});
end
  