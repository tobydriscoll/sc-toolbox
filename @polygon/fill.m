function H = fill(p,varargin)
%FILL   Plot a polygon with a filled interior.
%   FILL(P) plots the boundary of P in blue and fills the interior of the
%   polygon with gray. FILL(P,PROP1,VAL1,...) passes additional arguments
%   to the built-in FILL command.
%
%   See also FILL.

% Copyright 2003 by Toby Driscoll.
% $Id: fill.m 274 2004-05-27 13:11:21Z driscoll $

v = vertex(p);
vf = v(~isinf(v));
if any(isinf(v))
  v = vertex(truncate(p));
end

axlim = [min(real(vf)) max(real(vf)) min(imag(vf)) max(imag(vf))];
d = max([0.04,diff(axlim(1:2)),diff(axlim(3:4))]);
axlim(1:2) = mean(axlim(1:2)) + 0.54*[-1 1]*d;
axlim(3:4) = mean(axlim(3:4)) + 0.54*[-1 1]*d;

% Use defaults, but allow overrides and additional settings.
settings = { 0.75*[1 1 1],'edgecolor','b','linewidth',1.5, varargin{:} };
h = fill(real(v),imag(v),settings{:});

if ~ishold
  axis equal
  axis square
  axis(axlim)
end

if nargout > 0
  H = h;
end
