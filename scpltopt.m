function [nqpts,minlen,maxlen,maxrefn] = scpltopt(options)
%SCPLTOPT Parameters used by S-C plotting routines.
%   OPTIONS(1): Number of quadrature points per integration.
%               Approximately equals -log10(error).  Increase if plot
%               has false little zigzags in curves (default 4). 
%   OPTIONS(2): Minimum line segment length, as a proportion of the
%               axes box (default 0.005).
%   OPTIONS(3): Maximum line segment length, as a proportion of the
%               axes box (default 0.05).
%   OPTIONS(4): Max allowed number of adaptive refinements made to meet
%               other requirements (default 10).
%      
%   See also HPPLOT, DPLOT, DEPLOT, STPLOT, RPLOT, CRPLOT, CRRPLOT.

%   Copyright 1998 by Toby Driscoll.
%   $Id: scpltopt.m 298 2009-09-15 14:36:37Z driscoll $

user = options;
lenu = length(user);
options = zeros(1,4);
options(1:lenu) = user(1:lenu);
options = options + (options==0).*[5,.005,.02,16];

nqpts = options(1);
minlen = options(2);
maxlen = options(3);
maxrefn = options(4);
 