function f = composite(varargin)
%COMPOSITE Form a composition of maps.
%   F = COMPOSITE(F1,F2,...) defines F as the composite map obtained by
%   first applying F1, then F2, etc. Each member map may be an SCMAP,
%   SCMAPINV (inverse SC), INLINE function, or MOEBIUS map. Note that
%   composites not including an INLINE member may be inverted using INV.
%   
%   See also COMPOSITE/EVAL, COMPOSITE/INV, SCMAP, SCMAPINV, MOEBIUS

%   Copyright 2001 by Toby Driscoll.
%   $Id: composite.m 195 2002-09-10 19:11:44Z driscoll $

f.maps = {};

for n = 1:nargin
  map = varargin{n};
  switch class(map)
    case 'composite'
      m = members(map);
      f.maps = {f.maps{:}, m{:}};
    case {'moebius','diskmap','hplmap','extermap','stripmap','rectmap',...
	    'crdiskmap','crrectmap','riesurfmap','scmapinv'}
      f.maps{end+1} = map;
    case 'inline'
      if nargin(map) > 1
	error('Inline functions must have only one argument.')
      else
	f.maps{end+1} = map;
      end
    otherwise
      error('Object type ''%s'' not recognized.', class(map))
  end
end

f = class(f,'composite');

      
    