function [H,RE,IM] = stplot(w,beta,z,c,re,im,options)
%STPLOT Image of cartesian grid under Schwarz-Christoffel strip map.
%   STPLOT(W,BETA,Z,C) will adaptively plot the images under the
%   Schwarz-Christoffel exterior map of ten evenly spaced horizontal and
%   vertical lines in the upper half-plane. The abscissae of the
%   vertical lines will bracket the finite extremes of real(Z).  The
%   arguments are as in STPARAM.
%
%   STPLOT(W,BETA,Z,C,M,N) will plot images of M evenly spaced vertical
%   and N evenly spaced horizontal lines.  Horizontal lines are spaced
%   to bracket real(Z); vertical lines are evenly spaced between 0 and
%   1.
%
%   STPLOT(W,BETA,Z,C,RE,IM) will plot images of vertical lines whose
%   real parts are given in RE and horizontal lines whose imaginary
%   parts are given in IM.  Either argument may be empty.
%
%   STPLOT(W,BETA,Z,C,RE,IM,OPTIONS) allows customization of HPPLOT's
%   behavior.  See SCPLTOPT.
%
%   H = STPLOT(W,BETA,Z,C,...) returns a vector of handles to all the
%   curves drawn in the interior of the polygon.  [H,RE,IM] =
%   STPLOT(W,BETA,Z,C,...) also returns the abscissae and ordinates of
%   the lines comprising the grid.
%   
%   See also SCPLTOPT, STPARAM, STMAP, STDISP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: stplot.m 298 2009-09-15 14:36:37Z driscoll $

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
minre = min(real(z(~isinf(z))));
maxre = max(real(z(~isinf(z))));
if (length(re)==1) & (re == round(re))
  if re < 1
    re = [];
  elseif re < 2
    re = mean([minre,maxre]);
  else
    m = re;
    re = linspace(minre,maxre,m);
    dre = diff(re(1:2));
    re = linspace(minre-dre,maxre+dre,m);
  end
end
if (length(im)==1) & (im == round(im))
  if im < 1
    im = [];
  else
    m = im;
    im = linspace(0,1,m+2);
    im([1,m+2]) = [];
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
qdat = scqdata(beta,nqpts);
len = max(diff(get(ax(1),'xlim')),diff(get(ax(1),'ylim')));
minlen = len*minlen;
maxlen = len*maxlen;
axlim = axis;

color = 'k';

% Plot vertical lines...
linh = zeros(length(re),2);
for j = 1:length(re)
  % Start evenly spaced
  zp = re(j) + i*linspace(0,1,15).';
  new = logical(ones(size(zp)));
  wp = repmat(NaN,length(zp),1);

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
    neww = stmap(zp(new),w,beta,z,c,qdat);
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
	'xdata',re(j)*[1 1],'ydata',[0 1])
  end
  drawnow
end

% Plot horizontal lines...
x1 = min(-5,minre);
x2 = max(5,maxre);
wl = w(isinf(z) & z < 0);		% image of left end
wr = w(isinf(z) & z > 0);		% image of right end
linh1 = linh;
linh = zeros(length(im),2);
for j = 1:length(im)
  % Start evenly spaced
  zp = [-Inf linspace(x1,x2,15) Inf].' + i*im(j);
  new = logical(ones(size(zp)));
  new([1 end]) = 0;
  wp = repmat(NaN,length(zp),1);
  wp([1 end]) = [wl; wr];
  
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
    neww = stmap(zp(new),w,beta,z,c,qdat);
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
  
% Finish up

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

