function varargout = get(p,varargin)
%GET    Get polygon parameters.
%   ==> Deprecated. Please use ANGLE and VERTEX instead.
%
%   [VAL1,VAL2,...] = GET(P,'PROP1','PROP2',...) returns the values of the
%   polygon P corresponding to the requested properties. Valid properties
%   are: 
%   
%       vertex, angle

% Copyright 1999-2003 by Toby Driscoll.
% $Id: get.m 232 2003-01-09 14:52:16Z driscoll $

for j = 1:length(varargin)
  switch lower(varargin{j}(1:min(3,length(varargin{j}))))
   case 'ver'
    varargin{j} = p.vertex;
   case 'ang'
    varargin{j} = p.angle;
   otherwise
    warning(sprintf('Property ''%s'' not recognized.\n',varargin{j}))
    varargin{j} = [];
  end
end
