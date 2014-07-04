function [H,R2,THETA] = crplot(w,beta,cr,aff,wcfix,Q,R,theta,options)
%CRPLOT Image of polar grid under disk map in crossratio form.
%   CRPLOT(W,BETA,CR,AFF,WCFIX,Q) will adaptively plot the images under
%   the Schwarz-Christoffel crossratio disk map of ten evenly spaced
%   circles and rays in the unit disk.
%
%   CRPLOT(W,BETA,CR,AFF,WCFIX,Q,M,N) will plot images of M evenly
%   spaced circles and N evenly spaced rays.
%
%   CRPLOT(W,BETA,CR,AFF,WCFIX,Q,R,THETA) will plot images of circles
%   whose radii are given in R and rays whose arguments are given in
%   THETA. Either argument may be empty.
%
%   CRPLOT(W,BETA,CR,AFF,WCFIX,Q,R,THETA,OPTIONS) allows customization
%   of CRPLOT's behavior. See SCPLTOPT.
%
%   H = CRPLOT(...) returns a vector of handles to all the curves drawn
%   in the interior of the polygon. [H,R,THETA] = CRPLOT(...)  also
%   returns the moduli and arguments of the curves comprising the grid.
%       
%   See also SCPLTOPT, CRPARAM, CRFIXWC, CRMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crplot.m 226 2003-01-08 16:59:17Z driscoll $

% Parse input and initialize
w = w(:);
beta = beta(:);
n = length(w);

% Parse input
if nargin < 9
  options = [];
  if nargin < 8
    theta = [];
    if nargin < 7 
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
  R = linspace(0,1,m);
  R([1,m]) = [];
end
if (length(theta)==1) && (theta == round(theta))
  m = theta+1;
  theta = linspace(0,2*pi,m);
  theta(m) = [];
end

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

% Drawing parameters
[nqpts,minlen,maxlen,maxrefn] = scpltopt(options);
qdat = scqdata(beta,nqpts);
len = max(diff(get(ax(1),'xlim')),diff(get(ax(1),'ylim')));
minlen = len*minlen;
maxlen = len*maxlen;
axlim = axis;

color = 'k';

% Plot circles...
linh = gobjects(length(R),2);
for j = 1:length(R)
  % Start with evenly spaced theta
  tp = linspace(0,2*pi,20)';
  new = true(length(tp),1);
  wp = NaN*new;

  % The individual points will be shown as they are found
  linh(j,1) = animatedline('parent',ax(1),'color',color,'vis',vis{1},...
      'linestyle','none','marker','.','markersize',7);
  if draw2
    linh(j,2) = animatedline('parent',ax(2),'color',color,'vis',vis{2},...
	'linestyle','none','marker','.','markersize',7);
  end
  
  % Adaptively refine theta to get smooth curve
  iter = 0;
  while (any(new)) && (iter < maxrefn)
    drawnow
    neww = crmap(R(j)*exp(1i*tp(new)),w,beta,cr,aff,wcfix,Q,qdat);
    wp(new) = neww;
    iter = iter + 1;

    % Update the points to show progress
    addpoints(linh(j,1),real(wp(new)),imag(wp(new)))
    if draw2
      addpoints(linh(j,2),R(j)*cos(tp(new)),R(j)*sin(tp(new)))
    end
    drawnow update
    
    % Add points to zp where necessary to make smooth curves
    [tp,wp,new] = scpadapt(tp,wp,minlen,maxlen,axlim);
  end
  % Set the lines to be solid
  % To ensure smooth closure, we erase and redraw all the points
  % in order. 
  clearpoints(linh(j,1))
  addpoints(linh(j,1),real(wp),imag(wp));    
  set(linh(j,1),'marker','none','linestyle','-','user',R(j)*exp(1i*tp))
  if draw2
    % Replace the points with (hopefully) a smooth circle
    tp = linspace(0,2*pi,361);
    clearpoints(linh(j,2))
    addpoints(linh(j,2),R(j)*cos(tp),R(j)*sin(tp));
    set(linh(j,2),'marker','none','linestyle','-')
  end
  drawnow
end

% Plot radii...
linh1 = linh;
linh = gobjects(length(theta),2);
for j = 1:length(theta)
  % Start with evenly spaced radii
  Rp = linspace(0,1,12)';
  Rp = [Rp(1:11);.95;.98;Rp(12)];
  zp = Rp*exp(1i*theta(j));
  new = true(length(zp),1);
  wp = NaN*new;

  % The individual points will be shown as they are found
  linh(j,1) = animatedline('parent',ax(1),'color',color,'vis',vis{1},...
      'linestyle','none','marker','.','markersize',7);
  if draw2
    linh(j,2) = animatedline('parent',ax(2),'color',color,'vis',vis{2},...
	'linestyle','none','marker','.','markersize',7);
  end

  % Adaptively refine to make smooth curves
  iter = 0;
  while (any(new)) && (iter < maxrefn)
    drawnow
    neww = crmap(zp(new),w,beta,cr,aff,wcfix,Q,qdat);
    wp(new) = neww;
    iter = iter + 1;

    % Update the points to show progress
    addpoints(linh(j,1),real(wp(new)),imag(wp(new)))
    if draw2
        addpoints(linh(j,1),real(zp(new)),imag(zp(new)))
    end
    drawnow update
    
    % Add points to zp where necessary to make smooth curves
    [zp,wp,new] = scpadapt(zp,wp,minlen,maxlen,axlim);
  end
  % Set the lines to be solid
  clearpoints(linh(j,1))
  addpoints(linh(j,1),real(wp),imag(wp));    
  set(linh(j,1),'marker','none','linestyle','-','user',zp)
  if draw2
    % Replace the points with just the ends
    clearpoints(linh(j,2))
    addpoints(linh(j,2),[0 1]*cos(theta(j)),[0 1]*sin(theta(j)))
    set(linh(j,2),'marker','none','linestyle','-')
  end
  drawnow
end

linh = [linh1;linh];
if ~draw2
  linh = linh(:,1);
end

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
