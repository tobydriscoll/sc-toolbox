function w = eval(f,z)
%EVAL   Evaluate a composite map.
%   EVAL(F,Z) evaluates the composite map F at Z.

% Copyright 2001 by Toby Driscoll.
% $Id: eval.m 169 2001-07-20 15:19:13Z driscoll $

w = z;
for n = 1:length(f.maps)
  w = feval(f.maps{n},w);
end
