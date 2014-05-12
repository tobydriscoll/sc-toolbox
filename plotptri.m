function h = plotptri(w,Q,lab)
%PLOTPTRI Plot a polygon triangulation.
%   PLOTPTRI(W,E) plots the polygon W and the edges of a triangulation
%   whose edges are described by E. Instead of E, you may also pass the
%   quadrilateral graph structure, Q.
%       
%   If a nonempty third argument is given, the edges and vertices are
%   labeled by number.
%       
%   An output argument will be assigned a vector of handles to the edges
%   drawn.

%   Copyright 1998 by Toby Driscoll.
%   $Id: plotptri.m 298 2009-09-15 14:36:37Z driscoll $

% Parse input
if isstruct(Q)
  edge = Q.edge;
else
  edge = Q;
end

% Plot all edges
we = edge;
we(:) = w(edge);
han = plot(we,'--','color',[0 .5 0]);

turn_off_hold = ~ishold;
hold on

% Label if requested
if nargin > 2 & ~isempty(lab)
  [he,hl] = plotpoly(w,scangle(w),1);

  % Add circles at vertices
  %hll = findobj(hl,'type','line');
  %set(hll,'marker','o','markersize',get(gca,'defaultlinemarkersize'))

  % Edge labels
  for k=1:size(edge,2)
    mp = mean(w(edge(:,k)));
    text(real(mp),imag(mp),int2str(k),'color',[0,0,0],...
	'vert','mid','hor','cen')
  end
else
  % Unlabeled version
  plotpoly(w,scangle(w));
  plot(w,'o')
end

if turn_off_hold
  hold off
end

if nargout > 0
  h = han;
end
