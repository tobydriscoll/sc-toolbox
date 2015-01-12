function [H,R2,THETA] = deplot(w,beta,z,c,R,theta,options)
%DEPLOT Image of polar grid under Schwarz-Christoffel exterior map.
%   DEPLOT(W,BETA,Z,C) will adaptively plot the images under the
%   Schwarz-Christoffel exterior map of ten evenly spaced circles and
%   rays in the unit disk.  The arguments are as in DEPARAM.
%
%   DEPLOT(W,BETA,Z,C,M,N) will plot images of M evenly spaced circles
%   and N evenly spaced rays.
%
%   DEPLOT(W,BETA,Z,C,R,THETA) will plot images of circles whose radii
%   are given in R and rays whose arguments are given in THETA.  Either
%   argument may be empty.
%
%   DEPLOT(W,BETA,Z,C,R,THETA,OPTIONS) allows customization of DEPLOT's
%   behavior.  See SCPLTOPT.
%
%   H = DEPLOT(W,BETA,Z,C,...) returns a vector of handles to all the
%   curves drawn in the interior of the polygon.  [H,R,THETA] =
%   DEPLOT(W,BETA,Z,C,...) also returns the moduli and arguments of the
%   curves comprising the grid.
%   
%   See also SCPLTOPT, DEPARAM, DEMAP, DEDISP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: deplot.m 226 2003-01-08 16:59:17Z driscoll $

beta = beta(:);
z = z(:);
import sctool.*

% Parse input
if nargin < 7
  options = [];
  if nargin < 6
    theta = [];
    if nargin < 5 
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
if (length(R)==1) && (R == round(R))
  m = R+2;
  R = fliplr(linspace(.25,1,m));
  R([1,m]) = [];
end
if (length(theta)==1) && (theta == round(theta))
  m = theta+1;
  theta = linspace(0,2*pi,m);
  theta(m) = [];
end

autoscale = strcmp(get(gca,'xlimmode'),'auto') && ...
    strcmp(get(gca,'ylimmode'),'auto');
autoscale = autoscale | ~ishold;

% Prepare figure
fig = gcf;
figure(fig);

% If figure has two tagged axes, draw in both
ax = [findobj(fig,'tag','PhysicalAxes');...
      findobj(fig,'tag','CanonicalAxes')];
if length(ax)==2
  draw2 = true;
  vis = get(ax,'vis');
else
  draw2 = false;
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

axlim = axis;
if autoscale
  axlim(1:2) = axlim(1:2) + 0.25*diff(axlim(1:2))*[-1,1];
  axlim(3:4) = axlim(3:4) + 0.25*diff(axlim(3:4))*[-1,1];
  axis(axlim);
end
drawnow

% Enlarge clipping axes to avoid cutting off arcs
axlim(1:2) = axlim(1:2) + 0.25*diff(axlim(1:2))*[-1,1];
axlim(3:4) = axlim(3:4) + 0.25*diff(axlim(3:4))*[-1,1];

% Drawing parameters
[nqpts,minlen,maxlen,maxrefn] = scpltopt(options);
qdat = scqdata(beta,nqpts);
len = max(diff(get(ax(1),'xlim')),diff(get(ax(1),'ylim')));
minlen = len*minlen;
maxlen = len*maxlen;

color = 'k';

% Plot circles...
linh = gobjects(length(R),2);
for j = 1:length(R)
  % Start with evenly spaced theta
  tp = linspace(0,2*pi,20)';
  new = true(length(tp),1);
  wp = NaN(size(new));

  % The individual points will be shown as they are found
  if verLessThan('matlab','8.4')
      linh(j,1) = line(NaN,NaN,'parent',ax(1),'color',color,'vis',vis{1},...
          'linestyle','none','marker','.','markersize',7,'erasemode','none');
      if draw2
          linh(j,2) = line(NaN,NaN,'parent',ax(2),'color',color,'vis',vis{2},...
              'linestyle','none','marker','.','markersize',7,'erasemode','none');
      end
  else
      
      linh(j,1) = animatedline('parent',ax(1),'color',color,'vis',vis{1},...
          'linestyle','none','marker','.','markersize',7);
      if draw2
          linh(j,2) = animatedline('parent',ax(2),'color',color,'vis',vis{2},...
              'linestyle','none','marker','.','markersize',7);
      end     
  end

  % Adaptively refine theta to make smooth curve
  iter = 0;
  while (any(new)) && (iter < maxrefn)
    drawnow
    zp = R(j)*exp(1i*tp(new));
    neww = demap(zp,w,beta,z,c,qdat);
    wp(new) = neww;
    iter = iter + 1;
    
    % Update the points to show progress
    if verLessThan('matlab','8.4')
        set(linh(j,1),'xdata',real(wp),'ydata',imag(wp))
        if draw2
            set(linh(j,2),'xdata',R(j)*cos(tp),'ydata',R(j)*sin(tp))
        end
    else
        addpoints(linh(j,1),real(wp(new)),imag(wp(new)))
        if draw2
            addpoints(linh(j,2),R(j)*cos(tp(new)),R(j)*sin(tp(new)))
        end
    end
    drawnow

    % Add points to zp where necessary
    [tp,wp,new] = scpadapt(tp,wp,minlen,maxlen,axis);

  end
  % Set the lines to be solid
  if verLessThan('matlab','8.4')
      set(linh(j,1),'erasemode','back')
      set(linh(j,1),'marker','none','linestyle','-','user',R(j)*exp(1i*tp))
      if draw2
          % Replace the points with (hopefully) a smooth circle
          tp = linspace(0,2*pi,101);
          set(linh(j,2),'erasemode','back')
          set(linh(j,2),'marker','none','linestyle','-',...
              'xdata',R(j)*cos(tp),'ydata',R(j)*sin(tp))
      end
      
  else
      
      clearpoints(linh(j,1))
      addpoints(linh(j,1),real(wp),imag(wp));
      set(linh(j,1),'marker','none','linestyle','-','user',R(j)*exp(1i*tp))
      if draw2
          % Replace the points with (hopefully) a smooth circle
          tp = linspace(0,2*pi,361);
          clearpoints(linh(j,2))
          addpoints(linh(j,2),R(j)*cos(tp),R(j)*sin(tp))
          set(linh(j,2),'marker','none','linestyle','-')
      end
  end
  
  drawnow
  
end


% Plot radii...
linh1 = linh;
linh = gobjects(length(theta),2);
for j = 1:length(theta)
  Rp = [0 linspace(.2,1,14)]';
  zp = Rp*exp(1i*theta(j));
  new = true(length(zp),1);
  wp = NaN(size(new));
  new(1) = 0;
  wp(1) = Inf;

  % The individual points will be shown as they are found
  if verLessThan('matlab','8.4')
      linh(j,1) = line(NaN,NaN,'parent',ax(1),'color',color,'vis',vis{1},...
          'linestyle','none','marker','.','markersize',7,'erasemode','none');
      if draw2
          linh(j,2) = line(NaN,NaN,'parent',ax(2),'color',color,'vis',vis{2},...
              'linestyle','none','marker','.','markersize',7,'erasemode','none');
      end
  else
      linh(j,1) = animatedline('parent',ax(1),'color',color,'vis',vis{1},...
          'linestyle','none','marker','.','markersize',7);
      if draw2
          linh(j,2) = animatedline('parent',ax(2),'color',color,'vis',vis{2},...
              'linestyle','none','marker','.','markersize',7);
      end
  end

  % Adaptively refine to make smooth curve
  iter = 0;
  while (any(new)) && (iter < maxrefn)
    drawnow
    neww = demap(zp(new),w,beta,z,c,qdat);
    wp(new) = neww;
    iter = iter + 1;

    % Update the points to show progress
    if verLessThan('matlab','8.4')
        set(linh(j,1),'xdata',real(wp),'ydata',imag(wp))
        if draw2
            set(linh(j,2),'xdata',real(zp),'ydata',imag(zp))
        end
    else
        addpoints(linh(j,1),real(wp(new)),imag(wp(new)))
        if draw2
            addpoints(linh(j,2),real(zp(new)),imag(zp(new)))
        end
    end
    drawnow
 
    % Add points to zp where necessary
    [zp,wp,new] = scpadapt(zp,wp,minlen,maxlen,axis);
  end

  % Set the lines to be solid
  if verLessThan('matlab','8.4')
      set(linh(j,1),'erasemode','back')
      set(linh(j,1),'marker','none','linestyle','-','user',zp)
      if draw2
          % Replace the points with just the ends
          set(linh(j,2),'erasemode','back')
          set(linh(j,2),'marker','none','linestyle','-',...
              'xdata',[0 1]*cos(theta(j)),'ydata',[0 1]*sin(theta(j)))
      end
  else
      clearpoints(linh(j,1))
      addpoints(linh(j,1),real(wp),imag(wp));
      set(linh(j,1),'marker','none','linestyle','-','user',zp)
      if draw2
          % Replace the points with just the ends
          clearpoints(linh(j,2))
          addpoints(linh(j,2),[0 1]*cos(theta(j)),[0 1]*sin(theta(j)))
          set(linh(j,2),'marker','none','linestyle','-')
      end
  end
  drawnow
end

linh = [linh1;linh];
if ~draw2
  linh = linh(:,1);
end
if verLessThan('matlab','8.4'), set(linh,'erasemode','normal'),end

drawnow

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

