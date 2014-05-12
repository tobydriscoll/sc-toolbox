function ptsource(w,beta,z,c,ws,R,theta,options)
%PTSOURCE Field due to point source in a polygon.
%   PTSOURCE plots evenly spaced equipotential and force lines for a
%   point source located in a polygonal region.  This is equivalent to
%   the disk map with conformal center at the source.  With no arguments
%   the user draws the polygon and clicks the mouse at the source.
%
%   PTSOURCE(W,BETA) uses the polygon described by W and BETA.
%
%   PTSOURCE(W,BETA,Z,C) assmues that Z and C comprise the solution to
%   the disk mapping parameter problem, as returned by DPARAM.
%      
%   PTSOURCE(W,BETA,Z,C,WS) uses WS as the source location.
%
%   PTSOURCE(W,BETA,Z,C,WS,R,THETA,OPTIONS) uses the R, THETA, and
%   OPTIONS parameter as described in SCPLOTOPT.
%   
%   See also DFIXWC.

%   Copyright 1998 by Toby Driscoll.
%   $Id: ptsource.m 298 2009-09-15 14:36:37Z driscoll $

if nargin < 2
  [w,beta] = drawpoly;
end 
n = length(w);
if nargin < 8
  options = [];
  if nargin < 7
    theta = [];
    if nargin < 6
      R = [];
      if nargin < 5
	ws = [];
	if nargin < 4
	  z = [];
	end
      end
    end
  end
end

if isempty(z)
  [z,c] = dparam(w,beta);
end

if isempty(ws)
  plotpoly(w,beta)
  disp('Click mouse at source location.')
  [xc,yc] = ginput(1);
  ws = xc+i*yc;
end

[z,c] = dfixwc(w,beta,z,c,ws);
dplot(w,beta,z,c,R,theta,options);
