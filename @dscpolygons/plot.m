function plot(rgn,varargin)
%PLOT Plot a doubly connected region.
%   PLOT(RGN) plots the doubly connected region.
%   
%   PLOT(RGN,'num') or PLOT(RGN,'lab') also plots dots for the vertices and
%   numeric labels of inner polygon and outer polygon. For infinite vertices, 
%   two numeric labels are printed.

%   Written by Alfa Heryudono, 2003.

hold on;
if nargin >1
    if strcmp(varargin,'lab')==1 
        plot(rgn.p0,'lab');
        plot(rgn.p1,'lab');
    elseif strcmp(varargin,'num')==1 
        plot(rgn.p0,'num');
        plot(rgn.p1,'num');
    else
        plot(rgn.p0);
        plot(rgn.p1);
    end
else
    plot(rgn.p0);
    plot(rgn.p1);
end
hold off;
