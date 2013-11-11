function varargout = feval(varargin)
%FEVAL   Equivalent to EVAL.

if nargout
  varargout = cell(1,nargout);
  [varargout{:}] = eval(varargin{:});
else
  varargout{1} = eval(varargin{:});
end


