function options = scmapopt(varargin)
%SCMAPOPT Set options for SC maps.
%   options = SCMAPOPT('NAME1',VALUE1,'NAME2',VALUE2,...) creates an
%   options structure in which the named properties have the specified
%   values. Any unspecified properties have default values.  It is
%   sufficient to type only the leading characters that uniquely
%   identify the property; case is ignored for property names.
%   
%   SCMAPOPT(OLDOPTS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTS.
%   
%   SCMAPOPT(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS.  Any new properties
%   overwrite corresponding old properties.
%   
%   SCMAPOPT with no arguments displays all property names and their
%   possible values.
%   
%PROPERTIES:
% 
%Method: [ fsolve | {nesolve} ]
%   Which nonlinear equations solver to use. fsolve is used by default
%   if the Optimization Toolbox version 2 or higher is present;
%   otherwise, the nesolve suite packaged with this toolbox is
%   used. (NOT YET IMPLEMENTED--ALWAYS USES NESOLVE)
%   
%TraceSolution: [ full | {on} | off ]
%   How much progress information to show during and after the solution
%   to the parameter problem.
%   
%Tolerance: [ {1e-8} ]
%   Desired accuracy in the map. This may not be met exactly.
%   
%SolverMethod: [ linesearch | {trust} ]
%   Different strategies in the nonlinear solver that attempt to
%   globalize convergence. (Only affects nesolve.)
%   
%InitialGuess: [ [] ]
%   If nonempty, used as the initial guess for the nonlinear solver.


%
%WindowPopup: [ 'on' | {'off'} ]
%   The STRIPMAP, RECTMAP, and CRRECTMAP solvers may need additional
%   graphical input. If this option is 'on', this will be done in a new
%   figure window which is then destroyed.

%   Copyright 1998--2001 by Toby Driscoll.
%   $Id: scmapopt.m 298 2009-09-15 14:36:37Z driscoll $
%   Adpated by Toby Driscoll from ODESET:
%     Mark W. Reichelt and Lawrence F. Shampine, 5/6/94
%     Copyright (c) 1984-96 by The MathWorks, Inc.


% Print out possible values of properties.
if (nargin == 0) & (nargout == 0)
  fprintf('   TraceSolution: [ full | {on} | off ]\n');
  fprintf('       Tolerance: [ {1e-8} ]\n');
  fprintf('    SolverMethod: [ linesearch | {trustregion} ]\n');
  fprintf('    InitialGuess: [ [] ]\n');
%  fprintf('    WindowPopup:  [ on | {off} ]\n');
  fprintf('\n');
  return;
end

Names = {
    'TraceSolution'
    'Tolerance'
    'SolverMethod'
    'InitialGuess'
%    'WindowPopup'
    };
[m,n] = size(Names);
names = lower(char(Names));

% Defaults
Defaults = {
    'on'
    1e-8
    'trust'
    []
%    'off'
    };

% Combine all leading options structures o1, o2, ... in scmapopt(o1,o2,...).
options = cell2struct(Defaults,Names,1);
i = 1;
while i <= nargin
  arg = varargin{i};
  if isstr(arg)                         % arg is an option name
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument
    if ~isa(arg,'struct')
      error(sprintf(['Expected %s to be a string property name ' ...
	 'or\nan options structure created with SCMAPOPT.'], inputname(i)));
    end
    for j = 1:m
      val = getfield(arg,Names{j});
      if ~isequal(val,[])             % empty strings '' do overwrite
	options = setfield(options,Names{j},val);
      end
    end
  end
  i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
  error('Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
  arg = varargin{i};
    
  if ~expectval
    if ~isstr(arg)
      error(sprintf('Expected %s to be a string property name.', ...
	  inputname(i)));
    end
    
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error(sprintf('Unrecognized property name ''%s''.', arg));
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
        msg = sprintf('Ambiguous property name ''%s'' ', arg);
        msg = [msg '(' Names{j(1)}];
        for k = j(2:length(j))'
          msg = [msg ', ' Names{k}];
        end
        msg = sprintf('%s).', msg);
        error(msg);
      end
    end
    expectval = 1;                      % we expect a value next
    
  else
    options = setfield(options,Names{j},arg);
    expectval = 0;
      
  end
  i = i + 1;
end

if expectval
  error(sprintf('Expected value for property ''%s''.', arg));
end
