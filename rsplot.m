function [H,R2,THETA] = rsplot(w,beta,z,zb,c,R,theta,options)
%RSPLOT  Image of polar grid under Schwarz-Christoffel RS map.
%   RSPLOT(W,BETA,Z,ZB,C) will adaptively plot the images under the
%   Schwarz-Christoffel disk map of ten evenly spaced circles and rays
%   in the unit disk.  
%
%   RSPLOT(W,BETA,Z,ZB,C,M,N) will plot images of M evenly spaced circles
%   and N evenly spaced rays.
%
%   RSPLOT(W,BETA,Z,ZB,C,R,THETA) will plot images of circles whose radii
%   are given in R and rays whose arguments are given in THETA.  Either
%   argument may be empty.
%
%   RSPLOT(W,BETA,Z,ZB,C,R,THETA,OPTIONS) allows customization of DPLOT's
%   behavior.  See SCPLTOPT.
%
%   H = RSPLOT(W,BETA,Z,ZB,C,...) returns a vector of handles to all the
%   curves drawn in the interior of the polygon.  [H,R,THETA] =
%   RSPLOT(W,BETA,Z,ZB,C,...) also returns the moduli and arguments of the
%   curves comprising the grid.
%   
%   See also SCPLTOPT, RSPARAM, RSMAP.

%   Copyright 2002 by Toby Driscoll.
%   $Id: rsplot.m 298 2009-09-15 14:36:37Z driscoll $

w = w(:);
beta = beta(:);
n = length(w);
z = z(:);
zb = zb(:);

% Parse input
if nargin < 8
  options = [];
  if nargin < 7
    theta = [];
    if nargin < 6 
      R = [];
    end
  end
end

% Empty arguments default to 10
if isempty([R(:);theta(:)])
  R = 10;
  theta = 10;
end

% Integer arguments must be converted to specific values
if (length(R)==1) & (R == round(R))
  m = R+2;
  R = linspace(0,1,m);
  R([1,m]) = [];
end
if (length(theta)==1) & (theta == round(theta))
  m = theta+1;
  theta = linspace(0,2*pi,m);
  theta(m) = [];
end

% Prepare figure
fig = gcf;
figure(fig);

% If figure has two tagged axes, draw in both
ax = [findobj(fig,'tag','phydomain');...
      findobj(fig,'tag','candomain')];
if length(ax)==2
  draw2 = logical(1);
  vis = get(ax,'vis');
else
  draw2 = logical(0);
  ax = gca;
  vis = {'on'};
end

% Prepare axes
axes(ax(1))
turn_off_hold = ~ishold;
hp = plotpoly(w,beta);
set(hp,'vis',vis{1})
hold on
% For now, there is no need to draw the canonical domain. This is an
% awkward truce with the GUI.

% Drawing parameters
[nqpts,minlen,maxlen,maxrefn] = scpltopt(options);
qdat = scqdata(beta,nqpts);
len = max(diff(get(ax(1),'xlim')),diff(get(ax(1),'ylim')));
minlen = len*minlen;
maxlen = len*maxlen;
axlim = axis;

color = 'k';

% Plot circles...
linh = zeros(length(R),2);
for j = 1:length(R)
  % Start with evenly spaced theta
  tp = linspace(0,2*pi,20)';
  new = logical(ones(length(tp),1));
  wp = NaN*new;

  % The individual points will be shown as they are found
  linh(j,1) = line(NaN,NaN,'parent',ax(1),'color',color,'vis',vis{1},...
      'linestyle','none','marker','.','markersize',7,'erasemode','none');
  if draw2
    linh(j,2) = line(NaN,NaN,'parent',ax(2),'color',color,'vis',vis{2},...
	'linestyle','none','marker','.','markersize',7,'erasemode','none');
  end
  
  % Adaptively refine theta to make smooth curve
  iter = 0;
  while (any(new)) & (iter < maxrefn)
    drawnow
    zp = R(j)*exp(i*tp(new));
    neww = rsmap(zp,w,beta,z,zb,c,qdat);
    wp(new) = neww;
    iter = iter + 1;
    
    % Update the points to show progress
    set(linh(j,1),'xdata',real(wp),'ydata',imag(wp))
    if draw2
      set(linh(j,2),'xdata',R(j)*cos(tp),'ydata',R(j)*sin(tp))
    end

    % Add points to zp where necessary
    [tp,wp,new] = scpadapt(tp,wp,minlen,maxlen,axlim);

  end
  
  % Set the lines to be solid
  set(linh(j,1),'erasemode','back')
  set(linh(j,1),'marker','none','linestyle','-','user',R(j)*exp(i*tp))
  if draw2
    % Replace the points with (hopefully) a smooth circle
    tp = linspace(0,2*pi,101);
    set(linh(j,2),'erasemode','back')
    set(linh(j,2),'marker','none','linestyle','-',...
	'xdata',R(j)*cos(tp),'ydata',R(j)*sin(tp))
  end
  drawnow
end

% Plot radii...
linh1 = linh;
linh = zeros(length(theta),2);
for j = 1:length(theta)
  Rp = linspace(0,1,14)';
  zp = Rp*exp(i*theta(j));
  new = logical(ones(length(zp),1));
  wp = NaN*new;

  % The individual points will be shown as they are found
  linh(j,1) = line(NaN,NaN,'parent',ax(1),'color',color,'vis',vis{1},...
      'linestyle','none','marker','.','markersize',7,'erasemode','none');
  if draw2
    linh(j,2) = line(NaN,NaN,'parent',ax(2),'color',color,'vis',vis{2},...
	'linestyle','none','marker','.','markersize',7,'erasemode','none');
  end
 
  % Adaptively refine to make smooth curve
  iter = 0;
  while (any(new)) & (iter < maxrefn)
    drawnow
    neww = rsmap(zp(new),w,beta,z,zb,c,qdat);
    wp(new) = neww;
    iter = iter + 1;

    % Update the points to show progress
    set(linh(j,1),'xdata',real(wp),'ydata',imag(wp))
    if draw2
      set(linh(j,2),'xdata',real(zp),'ydata',imag(zp))
    end
    
    % Add points to zp where necessary
    [zp,wp,new] = scpadapt(zp,wp,minlen,maxlen,axlim);
  end

  % Set the lines to be solid
  set(linh(j,1),'erasemode','back')
  set(linh(j,1),'marker','none','linestyle','-','user',zp)
  if draw2
    % Replace the points with just the ends
    set(linh(j,2),'erasemode','back')
    set(linh(j,2),'marker','none','linestyle','-',...
	'xdata',[0 1]*cos(theta(j)),'ydata',[0 1]*sin(theta(j)))
  end
  drawnow
end

linh = [linh1;linh];
if ~draw2
  linh = linh(:,1);
end
set(linh,'erasemode','normal')

refresh

if turn_off_hold, hold off, end;
if nargout > 0
  H = linh;
  if nargout > 1
    R2 = R;
    if nargout > 2
      THETA = theta;
    end
  end
end 

