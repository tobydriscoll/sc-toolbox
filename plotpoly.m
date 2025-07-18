function [edgehan,lblhan] = plotpoly(w,beta,~)
%PLOTPOLY Plot a (generalized) polygon.
%   PLOTPOLY(W,BETA) plots the polygon whose vertices are in vector W
%   and whose turning angles are in BETA. Vertices at infinity are
%   permitted, but there must be at least two consecutive finite
%   vertices somewhere in W.
%
%   PLOTPOLY(W,BETA,1) also plots markers at vertices and adds numeric
%   labels.
%   
%   HE = PLOTPOLY(W,BETA) or [HE,HL] = PLOTPOLY(W,BETA,1) returns
%   handles to the edges and markers/labels.
%   
%   See also DRAWPOLY, MODPOLY.

%   Copyright 1998 by Toby Driscoll.
%   $Id: plotpoly.m 298 2009-09-15 14:36:37Z driscoll $

lw = 2.5*get(gcf,'defaultlinelinewidth');
n = length(w);
atinf = isinf(w);
wf = w(~atinf);

if nargin == 1
  beta = sctool.scangle(w);
end

% Decide whether to rescale axes
turn_off_hold = ~ishold;
autoscale = strcmp(get(gca,'xlimmode'),'auto') & ...
    strcmp(get(gca,'ylimmode'),'auto');
autoscale = autoscale | turn_off_hold;
if autoscale
  lim = [min(real(wf)),max(real(wf)),min(imag(wf)),max(imag(wf))];
  maxdiff = max(diff(lim(1:2)),diff(lim(3:4)));
  fac = .6 + .1*(any(atinf));
  lim(1:2) = mean(lim(1:2)) + fac*maxdiff*[-1,1];
  lim(3:4) = mean(lim(3:4)) + fac*maxdiff*[-1,1];
else
  lim = axis;
end
R = max(lim(2)-lim(1),lim(4)-lim(3));

% Renumber to start with two finite vertices
first = find(~atinf & ~atinf([2:n,1]), 1 );
if isempty(first) 
  error('There must be two consecutive finite vertices.')
end
renum = [first:n,1:first-1];
w = w(renum);
beta = beta(renum);
atinf = isinf(w);

% First edge
edgeh = gobjects(n,1);
lblh = gobjects( n, 2 );
edgeh(1) = plot(real(w(1:2)),imag(w(1:2)),'-','linewid',lw);
ang = angle(w(2)-w(1));
hold on
axis(lim)

% Remaining edges
j = 2;
while j <= n
  jp1 = rem(j,n)+1;

  % Draw marker/label 
  if nargin == 3
    theta = ang;
    % May need to modify position of label
    if any( abs([beta(j-1)-1 beta(jp1)-1]) < 3*eps)
      % Next to a slit, perturb label inside
      theta = theta + pi - (beta(j)+1)*pi/2;
    elseif abs(beta(j)) < 3*eps
      % For a "trivial" vertex, put number outside
      theta = theta - pi/2;
    end
    % Make label; markers will be added last
    pos = w(j) + .035*R*exp(1i*theta);
    %%lblh(j,1) = plot(real(w(j)),imag(w(j)),'.','markersize',12);
    lblh(j,2) = text(real(pos),imag(pos),int2str(renum(j)),...
	'ver','mid','hor','cen');
  end
  
  % Next edge(s)
  if ~atinf(jp1)
    % Bounded edge; straightforward
    edgeh(j) = plot(real(w([j jp1])),imag(w([j jp1])),'-','linewid',lw);
    ang = ang - pi*beta(j);
    j = j+1;
  else
    % Unbounded edge (first of two consecutive)
    ang = ang-pi*beta(j);
    z = [w(j);w(j)+R*exp(1i*ang)];
    edgeh(j) = plot(real(z),imag(z),'-','linewid',lw);
    % Make first label outside axes box
    if nargin == 3
      theta = ang;
      Rx = (lim(1:2) - real(w(j))) / (cos(theta)+eps*(cos(theta)==0));
      Ry = (lim(3:4) - imag(w(j))) / (sin(theta)+eps*(sin(theta)==0));
      RR = [Rx,Ry];
      pos = w(j) + (min(RR(RR>0))+.07*R)*exp(1i*theta);
      str = sprintf('%i (inf)',renum(j+1));
      lblh(j+1,1) = text(real(pos),imag(pos),str,'ver','mid','hor','cen');
    end

    % Second unbounded edge
    ang = ang-pi*beta(jp1);
    z = [w(rem(j+1,n)+1)-R*exp(1i*ang);w(rem(j+1,n)+1)];
    edgeh(j+1) = plot(real(z),imag(z),'-','linewid',lw);
    if nargin == 3
      theta = ang + pi;
      Rx = (lim(1:2) - real(z(2))) / (cos(theta)+eps*(cos(theta)==0));
      Ry = (lim(3:4) - imag(z(2))) / (sin(theta)+eps*(sin(theta)==0));
      RR = [Rx,Ry];
      pos = z(2) + (min(RR(RR>0))+.07*R)*exp(1i*theta);
      str = sprintf('%i (inf)',renum(j+1));
      lblh(j+1,2) = text(real(pos),imag(pos),str,'ver','mid','hor','cen');
    end

    % We've done two
    j = j+2;
  end
end

% Last item: label for first point
if nargin == 3
  theta = ang;
  if any( abs([beta(n)-1 beta(1) beta(2)-1]) < 3*eps)
    theta = theta + pi -(beta(1)+1)*pi/2;
  end
  pos = w(1) + .035*R*exp(1i*theta);
  %%lblh(1,1) = plot(real(w(1)),imag(w(1)),'.','markersize',12);
  lblh(1,2) = text(real(pos),imag(pos),int2str(renum(1)),...
      'ver','mid','hor','cen');
end

% Plot markers last, to keep them "on top"
if nargin==3
  for j = find(~atinf(:)')
    lblh(j,1) = plot(real(w(j)),imag(w(j)),'.','markersize',12);
  end
  set(lblh,'color',get(edgeh(1),'color'))
end

% Clean up
set(gca,'dataaspectratio',[1 1 1])
set(gca,'plotboxaspectratio',[1 1 1])

if turn_off_hold
  hold off
end 

% Flag for the GUI to grab onto.
set(edgeh,'tag','PolygonSide');

% Output args
if nargout 
  edgehan(renum) = edgeh;
  if nargout > 1
    lblhan = zeros(n,2);
    lblhan(renum,:) = lblh;
  end
end

