function varargout =  parseopt(options)

%   Copyright 1998--2001 by Toby Driscoll.
%   $Id: parseopt.m 199 2002-09-13 18:54:27Z driscoll $

% There are 2 allowable inputs: old-fashioned array and newfangled
% structure. Handle both here.

if ~isstruct(options)
  user = options;
  lenu = length(user);
  options = zeros(1,3);
  options(1:lenu) = user(1:lenu);
  options = options + (options==0).*[0,1e-8,2];
  
  trace = options(1);
  tol = options(2);
  method = options(3);
  varargout = { trace,tol,method };
else
  switch(options.TraceSolution)
   case 'full'
    trace = 2;
   case 'on'
    trace = 1;
   otherwise
    trace = 0;
  end
  tol = options.Tolerance;
  method = strmatch(options.SolverMethod,{'line','trust'});
  %%newwindow = strcmp(options.WindowPopup,'on');
  varargout = { trace, tol, method };
end
