function fi = inv(f)
%INV    Invert a composite map if possible.
%   INV(F) will return a composite that is the inverse of F. However,
%   composites using INLINE maps cannot be inverted.

%   Copyright 2001 by Toby Driscoll.
%   $Id: inv.m 170 2001-07-20 15:19:52Z driscoll $

N = length(f.maps);
list = cell(1,N);
for n = 1:N
  m = N+1-n;
  if ~isa(f.maps{m},'inline')
    list{n} = inv(f.maps{m});
  else
    error('Can''t invert INLINE maps.')
  end
end
fi = composite(list{:});

