function [H,RE,IM] = hpplot(w,beta,z,c,re,im,options)
%HPPLOT Image of cartesian grid under Schwarz-Christoffel half-plane map.
%   HPPLOT(W,BETA,Z,C) will adaptively plot the images under the
%   Schwarz-Christoffel exterior map of ten evenly spaced horizontal and
%   vertical lines in the upper half-plane. The abscissae of the
%   vertical lines will bracket the finite extremes of Z.  The arguments
%   are as in HPPARAM.
%
%   HPPLOT(W,BETA,Z,C,M,N) will plot images of M evenly spaced vertical
%   and N evenly spaced horizontal lines.  The spacing will be the same
%   in both directions.
%
%   HPPLOT(W,BETA,Z,C,RE,IM) will plot images of vertical lines whose
%   real parts are given in RE and horizontal lines whose imaginary
%   parts are given in IM.  Either argument may be empty.
%
%   HPPLOT(W,BETA,Z,C,RE,IM,OPTIONS) allows customization of HPPLOT's
%   behavior.  See SCPLTOPT.
%
%   H = HPPLOT(W,BETA,Z,C,...) returns a vector of handles to all the
%   curves drawn in the interior of the polygon.  [H,RE,IM] =
%   HPPLOT(W,BETA,Z,C,...) also returns the abscissae and ordinates of
%   the lines comprising the grid.
%   
%   See also SCPLTOPT, HPPARAM, HPMAP, HPDISP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: hpplot.m 298 2009-09-15 14:36:37Z driscoll $

n = length(w);
w = w(:);
beta = beta(:);
z = z(:);

% Parse input
if nargin < 7
  options = [];
  if nargin < 6
    im = [];
    if nargin < 5
      re = [];
    end
  end
end

% Empty arguments default to 10
if isempty([re(:);im(:)])
  re = 10;
  im = 10;
end

% Integer arguments must be converted to specific values
if (length(re)==1) & (re == round(re))
  if re < 1
    re = [];
  elseif re < 2
    re = mean(z([1,n-1]));
  else
    m = re;
    re = linspace(z(1),z(n-1),m);
    dre = diff(re(1:2));
    re = linspace(z(1)-dre,z(n-1)+dre,m);
  end
end
if (length(im)==1) & (im == round(im))
  if length(re) < 2
    im = linspace(0,4,im+1);
    im(1) = [];
  else
    im = mean(diff(re))*(1:im);
  end
end

% Prepare figure
fig = gcf;
figure(fig);

% If figure has two tagged axes, draw in both
ax = [findobj(fig,'tag','PhysicalAxes');...
      findobj(fig,'tag','CanonicalAxes')];
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
len = max(diff(get(ax(1),'xlim')),diff(get(ax(1),'ylim')));
minlen = len*minlen;
maxlen = len*maxlen;
if any(isinf(z))
  qdat = scqdata(beta(1:n-1),nqpts);
else
  qdat = scqdata(beta,nqpts);
end
axlim = axis;

color = 'k';

% Plot vertical lines...
y2 = max(z(n-1),10);
linh = gobjects(length(re),2);
for j = 1:length(re)
  % Start evenly spaced
  zp = re(j) + i*[linspace(0,y2,15) Inf].';
  new = logical(ones(size(zp)));
  new(end) = 0;
  wp = repmat(NaN,length(zp),1);
  wp(end) = w(n);

  % The individual points will be shown as they are found
  linh(j,1) = line(NaN,NaN,'parent',ax(1),'color',color,'vis',vis{1},...
      'linestyle','none','marker','.','markersize',7,'erasemode','none');
  if draw2
    linh(j,2) = line(NaN,NaN,'parent',ax(2),'color',color,'vis',vis{2},...
	'linestyle','none','marker','.','markersize',7,'erasemode','none');
  end
  
  % Adaptive refinement to make smooth curve
  iter = 0;
  while (any(new)) & (iter < maxrefn)
    drawnow
    neww = hpmap(zp(new),w,beta,z,c,qdat);
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
    % Replace the points with the endpoints
    set(linh(j,2),'erasemode','back')
    set(linh(j,2),'marker','none','linestyle','-',...
	'xdata',re(j)*[1 1],'ydata',[0 imag(zp(end-1))])
  end
  drawnow
end

z1 = min(-10,z(1));
z2 = max(40,z(n-1));
linh1 = linh;
linh = gobjects(length(im),2);
for j = 1:length(im)
  % Start evenly spaced
  zp = [-Inf linspace(z1,z2,15) Inf].' + i*im(j);
  new = logical(ones(size(zp)));
  new([1 end]) = 0;
  wp = repmat(NaN,length(zp),1);
  wp([1 end]) = w(n);
  
  % The individual points will be shown as they are found
  linh(j,1) = line(NaN,NaN,'parent',ax(1),'color',color,'vis',vis{1},...
      'linestyle','none','marker','.','markersize',7,'erasemode','none');
  if draw2
    linh(j,2) = line(NaN,NaN,'parent',ax(2),'color',color,'vis',vis{2},...
	'linestyle','none','marker','.','markersize',7,'erasemode','none');
  end
 
  % Adaptive refinement to make smooth curve
  iter = 0;
  while (any(new)) & (iter < maxrefn)
    drawnow
    neww = hpmap(zp(new),w,beta,z,c,qdat);
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
    % Replace the points with the endpoints
    set(linh(j,2),'erasemode','back')
    set(linh(j,2),'marker','none','linestyle','-',...
	'xdata',real(zp([2 end-1])),'ydata',im(j)*[1 1]);
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
    RE = re;
    if nargout > 2
      IM = im;
    end
  end
end 

