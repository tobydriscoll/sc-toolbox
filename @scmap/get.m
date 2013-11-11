function varargout = get(map,varargin)
%GET    Get map parameters.
%   Each SC map stores data needed to compute with the map. GET is a
%   generic access to these data. 
%
%   [VAL1,VAL2,...] = GET(F,'PROP1','PROP2',...) returns the values of
%   the map F corresponding to the requested properties. Valid properties
%   vary by map type. Use GET(F) to see a list for the type associated
%   with F. Note that field name abbreviations are NOT allowed.

% Copyright 2003 by Toby Driscoll.
% $Id: get.m 239 2003-01-15 16:02:27Z driscoll $

if nargin==1   % typeout available fields

  % Get names.
  names = fieldnames(map);
  % Strip out the generic scmap (this class!).
  names( strmatch('scmap',names) ) = [];
  % Add in the always-present 'polygon' and 'options'.
  names = { 'polygon', 'options', names{:} };

  fprintf('\nAvailable fields in class %s:\n',class(map))
  fprintf('  %s\n',names{:})
  fprintf('\n')

else   % get values

  for j = 1:length(varargin)
    switch lower(varargin{j}(1:min(3,length(varargin{j}))))
     case 'pol'
      varargout{j} = polygon(map);
     case 'opt'
      varargout{j} = options(map);
     otherwise
      try
        % This is an affront to OO programming! However, it's very
        % convenient. 
        m = struct(map);
        varargout{j} = m.(lower(varargin{j}));
      catch
        warning(sprintf('Field ''%s'' not recognized.\n',varargin{j}))
        varargout{j} = [];
      end
    end
  end
  
end
