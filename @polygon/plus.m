function r = plus(p,q)
%   Translate a polygon, or add the vertices of two polygons.

%   Copyright 1998-2003 by Toby Driscoll.
%   $Id: plus.m 245 2003-03-03 16:24:51Z driscoll $

if isa(q,'polygon')
  tmp = p;
  p = q;
  q = tmp;
end

switch(class(q))
 case 'polygon'
  if length(q)~=length(p)
    error('Polygons must have the same length to be added.')
  elseif isinf(p) || isinf(q)
    error('Only finite polygons may be added.')
  end
  r = polygon( vertex(p) + vertex(q) );
 case 'double'
  if length(q) > 1 && length(q)~=length(p)
    error(['Only a scalar or identical-length vector may be added to a ' ...
           'polygon.'])
  end
  r = polygon( vertex(p) + q(:) );
end
