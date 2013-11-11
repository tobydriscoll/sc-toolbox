function v = parameters(M)
%PARAMETERS Return a structure of the Schwarz-Christoffel map parameters.

%   Copyright 2002 by Toby Driscoll.
%   $Id: parameters.m 201 2002-09-13 19:46:08Z driscoll $

v.branch = M.branch;
v.prevertex = M.prevertex;
v.prebranch = M.prebranch;
v.constant = M.constant;
