function [w,beta,indx] = modpoly(w,beta)
%MODPOLY Modify a polygon.
%   [WNEW,BETANEW] = MODPOLY(W,BETA) plots the polygon given by W and
%   BETA and allows the user to change it with the mouse.  At the start,
%   MODPOLY allows you to move vertices.  Move the cursor over a vertex
%   you want to move, hold down the left mouse button, drag it to its
%   new location, and release.  The vertex changes color and affected
%   sides become dashed as you move the vertex.
%   
%   To delete a vertex, press the Delete button.  The pointer will
%   change to a fleur.  After your next click and release, the selected
%   vertex will be deleted and you will return to movement mode.  To
%   cancel a requested deletion, press Delete again.
%   
%   The Add button works similarly.  To add, press the button and then
%   click and release on a polygon side.  A vertex will be added to the
%   middle of the side, the polygon is redrawn, and you return to
%   movement mode.
%   
%   Infinite vertices cannot be moved(!), deleted, or added.  When
%   moving the neighbor of an infinite vertex, the angle at infinity is
%   kept constant.  When you delete a neighbor of infinity, the turn at
%   the deleted vertex is lost and the angle at infinity changes.  You
%   cannot delete a vertex with two infinite neighbors.  When you add a
%   vertex to an infinite side, the new vertex appears at a "reasonable"
%   distance from its finite neighbor.
%   
%   [WNEW,BETANEW,IDX] = MODPOLY(W,BETA) also returns an index vector to
%   help keep track of additions and deletions.  IDX has the same length
%   as WNEW, and if IDX(J) is an integer, it gives the index that
%   WNEW(J) had in the original W.  If WNEW(J) was added, then IDX(J) is
%   NaN.
%   
%   Note: MODPOLY makes no attempt to keep the polygon "legal."  You can
%   easily create things which are not polygons, or change infinite
%   vertices into unrecognized finite ones.
%   
%   See also DRAWPOLY, PLOTPOLY.

%   Copyright 1998 by Toby Driscoll.
%   $Id: modpoly.m 298 2009-09-15 14:36:37Z driscoll $

global sc_hs sc_hv sc_k sc_w sc_beta sc_idx sc_color
n = length(sc_w); 

if nargin == 1 & isstr(w)
  % Pass through to callback
  cback(w);
  return
end

% Set up figure, axes, etc.
fig = gcf;
figure(fig);

% Save current properties to be restored
figprops = {'unit';'pos';'resizefcn';'pointer';'windowbuttondown';...
    'windowbuttonup';'pointer';'interruptible'};
figstate = get(fig,figprops);
axprops = {'unit';'pos';'xgrid';'ygrid';'interruptible'};
axstate = get(gca,axprops);

% Get figure position and screen size
set(fig,'unit','point')
figpos = get(fig,'pos');
unit = get(0,'unit');
set(0,'unit','point')
maxpos = get(0,'screensize');
set(0,'unit',unit)

% Draw polygon and initialize global vars
turn_off_hold = ~ishold;
n = length(w);
sc_w = w;				% global vertices
sc_beta = beta;			% global angles
cla
sc_hs = plotpoly(w,beta);			% side handles
% One side of a crack should be a different color.
crack = find(beta==1);
set(sc_hs(crack),'color',get(sc_hs(1),'color')/2)
hold on
sc_hv = repmat(NaN,n,1);			% vertex handles
for j=find(~isinf(w))'
  sc_hv(j) = plot(real(w(j)),imag(w(j)),'.','markersize',20);
end
set(sc_hv(rem(crack,n)+1),'color',get(sc_hv(1),'color')/2)
drawnow
axlim = axis;
sc_idx = (1:length(w))';		% indices
set(gcf,'pointer','circle')
set(gcf,'interruptible','off','busyaction','cancel')
set(gca,'interruptible','off','busyaction','cancel')

% Set up nice initial grid spacing
space = diff(axlim(1:2))/16;
space = min(0.5, ceil(space*2)/2);
axlim = round(axlim*2)/2;
axis normal
axis(axlim)

set(sc_hs,'erasemode','xor','buttondown','modpoly(''down'');')
set(sc_hv(ishandle(sc_hv)),'erasemode','xor','buttondown','modpoly(''down'');')
drawnow 

% Create uicontrols
[frame,control] = mp_make_widgets(fig);
mp_resize_widgets(fig);
set(fig,'resizefcn','modpoly(''resize'');')

%set(gcf,'windowbuttondown','modpoly(''down'');')

% Run in place until finished
uiwait(fig)

set(gcf,'windowbuttondown','')

% Recover new info and clean up
if get(findobj(gcf,'tag','mp_done'),'user') > 0
  w = sc_w;
  beta = sc_beta;
end
if turn_off_hold
  hold off
end
delete(control)
delete(frame)
delete(sc_hs)
delete(sc_hv(ishandle(sc_hv)))
set(gcf,figprops,figstate)
set(gca,axprops,axstate)

axis auto
plotpoly(w,beta)
indx = sc_idx;
clear global sc_hs sc_hv sc_k sc_w sc_beta sc_idx


function cback(action)

ptr = get(gcf,'pointer');
global sc_hs sc_hv sc_k sc_w sc_beta sc_idx
n = length(sc_w);
  
switch action
    
case 'down'				% button down
  h = get(gcf,'currentobj');
  % Act only if h is a line object
  if strcmp(get(h,'type'),'line')
    if ~strcmp(ptr,'crosshair')		% move or delete
      if strcmp(get(h,'marker'),'.')	% vertex?
	sc_k = find(h==sc_hv);
	colr = get(gca,'colororder');
	%set(h,'color',colr(2,:))
	set(sc_hs([sc_k,rem(sc_k-2+n,n)+1]),'linesty','--')
	if strcmp(ptr,'circle')
	  % Mouse movement needed only when moving vertices
	  set(gcf,'windowbuttonmotion','modpoly(''move'');')
	end
	set(gcf,'windowbuttonup','modpoly(''up'');')
      end
    else				% insert
      if strcmp(get(h,'linesty'),'-')	% edge?
	sc_k = find(h==sc_hs);
	set(sc_hs(sc_k),'linesty','--')
	set(gcf,'windowbuttonup','modpoly(''up'');')
      end	
    end
  end
  
  
case 'move'				% mouse move
  z = get(gca,'currentpoint');
  % Snap to grid if requested.
  if get(findobj(gcf,'tag','snap'),'value');
    grid = str2num(get(findobj(gcf,'tag','grid_space'),'string'))*[1 1];
    axlim = axis;
    minxy = axlim([1 3]);
    z(1,1:2) = minxy + grid.*(round((z(1,1:2)-minxy)./grid));
  end
  
  k = sc_k;
  set(sc_hv(k),'xd',z(1,1),'yd',z(1,2))
  
  % Must handle case of infinite predecessor/successor separately.
  j = rem(k,n)+1;			% successor
  if isinf(sc_w(j))
    xd = get(sc_hs(k),'xd');
    yd = get(sc_hs(k),'yd');
    phi = atan2(diff(yd),diff(xd));
    r = sqrt(diff(xd)^2+diff(yd)^2);
    y = sc_w(k) + [0,r*exp(i*phi)];
  else
    y = [z(1,1)+i*z(1,2),sc_w(j)];
    phi = angle(-diff(y)/diff(sc_w([j,k])));
    sc_beta(k) = sc_beta(k)-phi/pi;
    sc_beta(j) = sc_beta(j)+phi/pi;
  end
  set(sc_hs(k),'xd',real(y),'yd',imag(y))
  j = rem(k-2+n,n)+1;			% predecessor
  if isinf(sc_w(j))
    xd = get(sc_hs(j),'xd');
    yd = get(sc_hs(j),'yd');
    phi = atan2(-diff(yd),-diff(xd));
    r = sqrt(diff(xd)^2+diff(yd)^2);
    y = sc_w(k) + [r*exp(i*phi),0];
  else
    y = [sc_w(j),z(1,1)+i*z(1,2)];
    phi = angle(diff(y)/diff(sc_w([j,k])));
    sc_beta(k) = sc_beta(k)+phi/pi;
    sc_beta(j) = sc_beta(j)-phi/pi;
  end
  % Make change effective
  set(sc_hs(j),'xd',real(y),'yd',imag(y))
  drawnow
  sc_w(k) = z(1,1)+i*z(1,2);
  
case 'up'				% button up
  set(sc_hs([sc_k,rem(sc_k-2+n,n)+1]),'linesty','-')
  set(gcf,'windowbuttonup','')
  if strcmp(ptr,'circle')
    % Moved a vertex.  Just clean up.
    colr = get(gca,'colororder');
    set(sc_hv(sc_k),'color',colr(1,:))
    set(gcf,'windowbuttonmotion','')
  elseif strcmp(ptr,'fleur') 		% Delete...
	colr = get(gca,'colororder');
    set(sc_hv(sc_k),'color',colr(1,:))
    set(gcf,'pointer','circle')
    idx = rem(sc_k+(-2:1)+n-1,n)+1;	% 2 back, here, and 1 forward
    infb = isinf(sc_w(idx(2)));
    infa = isinf(sc_w(idx(4)));
    if n <= 3 | (infa & infb)
      return				% do nothing
    elseif ~infb & ~infa 
      % Finite neighborhs; easy.
      v = get(sc_hs(idx(1)),'xdata')+i*get(sc_hs(idx(1)),'ydata');
      v(3:4) = get(sc_hs(idx(4)),'xdata')+i*get(sc_hs(idx(4)),'ydata');
      b = scangle(v);
      sc_beta(idx([2,4])) = b(2:3);
      v = v(2:3);
    else
      % An infinite neighbor
      axlim = axis;
      r = sqrt(diff(axlim(1:2))^2+diff(axlim(3:4))^2);
      x = get(sc_hs(sc_k*infb + idx(2)*infa),'xdata');
      y = get(sc_hs(sc_k*infb + idx(2)*infa),'ydata');
      ang = atan2(diff(y),diff(x)) + (pi*infb);
      j = (idx(2)*infb) + (idx(4)*infa);
      sc_beta(j) = sc_beta(j) + sc_beta(sc_k);
      if infb
	v = sc_w(idx(4)) + [1.1*r*exp(i*ang),0];
      else
	v = sc_w(idx(2)) + [0,1.1*r*exp(i*ang)];
      end
    end
    
    set(sc_hs(idx(2)),...
	'xdata',real(v),'ydata',imag(v),'linesty','-')
    delete(sc_hv(sc_k))
    sc_hv(sc_k) = [];
    sc_w(sc_k) = [];
    delete(sc_hs(sc_k))
    sc_hs(sc_k) = [];
    sc_beta(sc_k) = [];
    sc_idx(sc_k) = [];
    
  elseif strcmp(ptr,'crosshair') 	% Add
    [wn,bn] = scaddvtx(sc_w,sc_beta,sc_k);
    sc_w = wn(:);
    sc_beta = bn(:);
    cla
    axlim = axis;
    sc_hs = plotpoly(wn,bn);
    axis(axlim)
    axis normal
    cback('snap');
    sc_hv = repmat(NaN,n+1,1);
    for j=find(~isinf(wn))'
      sc_hv(j) = plot(real(wn(j)),imag(wn(j)),'.','markersize',20);
    end
    set(sc_hs,'erasemode','xor','buttondown','modpoly(''down'');')
    set(sc_hv(ishandle(sc_hv)),'erasemode','xor',...
	'buttondown','modpoly(''down'');')
    drawnow 
    sc_idx = [sc_idx(1:sc_k);NaN;sc_idx(sc_k+1:n)];
    set(gcf,'pointer','circle')
    
  end
  
  
case 'delete'		 		% toggle delete state
  if ~strcmp(ptr,'fleur')
    set(gcf,'pointer','fleur')
  else
    set(gcf,'pointer','circle')
  end
  
  
case 'add'				% toggle add state
  if ~strcmp(ptr,'crosshair')
    set(gcf,'pointer','crosshair')
  else
    set(gcf,'pointer','circle')
  end
  
  
case 'done'
  set(findobj(gcf,'tag','mp_done'),'user',1)
  drawnow
  uiresume(gcf)
  
case 'cancel'
  set(findobj(gcf,'tag','mp_done'),'user',-1)
  drawnow
  uiresume(gcf)
  
case 'snap'
  fig = gcf;
  snap = findobj(fig,'tag','snap');
  if get(snap,'value')
    space = str2num(get(findobj(fig,'tag','grid_space'),'string'));
    axlim = axis;
    x = axlim(1):space:axlim(2);
    y = axlim(3):space:axlim(4);
    set(gca,'xticklabelmode','auto','yticklabelmode','auto')
    set(gca,'xtick',x,'ytick',y)
    % For clarity, keep only about eight of the labels.
    N = length(x);
    keep = [1,3:ceil((N-1)/8):N-2,N];
    xl = get(gca,'xticklabel');
    p = min(size(xl,2),4);
    xlnew = setstr(ones(N+1,1)*blanks(p));
    xlnew(keep,:) = xl(keep,1:p);
    N = length(y);
    keep = [1,3:ceil((N-1)/8):N-2,N];
    yl = get(gca,'yticklabel');
    p = min(size(yl,2),4);
    ylnew = setstr(ones(N+1,1)*blanks(p));
    ylnew(keep,:) = yl(keep,1:p);
    % Make it so.
    set(gca,'xticklabel',xlnew,'yticklabel',ylnew)
    set(gca,'xgrid','on','ygrid','on')
  else
    set(gca,'xtickmode','auto','xticklabelmode','auto','xgrid','off')
    set(gca,'ytickmode','auto','yticklabelmode','auto','ygrid','off')
  end    
  drawnow    
  
  
case 'xlim'
  xlim = str2num(get(findobj(gcf,'tag','xlimits'),'string'));
  set(gca,'xlim',xlim)
  cback('snap');
  
  
case 'ylim'
  ylim = str2num(get(findobj(gcf,'tag','ylimits'),'string'));
  set(gca,'ylim',ylim)
  cback('snap');
  
case 'resize'
  mp_resize_widgets(gcbo);
  
end


function [frame,widget] = mp_make_widgets(fig)
bgcolor = .7*[1 1 1];
set(fig,'defaultuicontrolbackground',bgcolor,'defaultuicontrolunits','point')

% Create uicontrol frames
frame(1) = uicontrol('style','frame','tag','mp_settings');
frame(2) = uicontrol('style','frame','tag','mp_actions');

axlim = axis;
space = 0.5;

widget(1) = uicontrol('style','check',...
    'string','Snap to grid',...
    'call','modpoly(''snap'');',...
    'hor','c',...
    'tag','snap');
widget(2) = uicontrol('style','text',...
    'string','Grid spacing:',...
    'hor','l');
widget(3) = uicontrol('style','edit',...
    'string',sprintf('%.3g',space),...
    'call','modpoly(''snap'');',...
    'tag','grid_space');

widget(4) = uicontrol('style','text',...
    'string','x limits:',...
    'hor','r');
widget(5) = uicontrol('style','edit',...
    'string',sprintf('[%.3g %.3g]',axlim(1:2)),...
    'hor','l',...
    'tag','xlimits',...
    'call','modpoly(''xlim'');');
widget(6) = uicontrol('style','text',...
    'string','y limits:',...
    'hor','r');
widget(7) = uicontrol('style','edit',...
    'string',sprintf('[%.3g %.3g]',axlim(3:4)),...
    'hor','l',...
    'tag','ylimits',...
    'call','modpoly(''ylim'');');

widget(8) = uicontrol('style','push',...
    'string','Add',...
    'call','modpoly(''add'');');
widget(9) = uicontrol('style','push',...
    'string','Delete',...
    'call','modpoly(''delete'');');
widget(10) = uicontrol('style','push',...
    'string','Done',...
    'call','modpoly(''done'');',...
    'tag','mp_done');
widget(11) = uicontrol('style','push',...
    'string','Cancel',...
    'call','modpoly(''cancel'');');

set(widget(10),'user',0)
set(frame(1),'userdata',widget(1:7))
set(frame(2),'userdata',widget(8:11))

function mp_resize_widgets(fig)

figpos = get(fig,'pos');

% Frames
f(1) = findobj(fig,'tag','mp_settings');
f(2) = findobj(fig,'tag','mp_actions');

% Minimum sizes
figpos(3) = max(figpos(3),376);
figpos(4) = max(figpos(4),170);

set(fig,'pos',figpos)
set(f,{'pos'},{[0 0 figpos(3) 72];[figpos(3)-80 72 80 figpos(4)-72]})

% Axes
% Axes area needs a border of at least 24 pt on each side
usable = figpos(3:4) - [80 72];
size = min(usable - 48);
margin = (usable - size)/2 + [0 72];
set(gca,'unit','point','pos',[margin(1:2) size size])

% Settings
h = get(f(1),'userdata');
pos = {
[4 48 100 20];...
[4 28 60 14];...
[62 26 40 18];...
[120 50 48 14];...
[172 49 54 18];...
[120 28 48 14];...
[172 27 54 18]
};
set(h,{'pos'},pos)

% Actions
offset = [figpos(3)-76 figpos(4)-94 0 0];
h = get(f(2),'userdata');
pos = {
[0 72 72 18]+offset;...
[0 52 72 18]+offset;...
[0 20 72 18]+offset;...
[0  0 72 18]+offset
};
set(h,{'pos'},pos)

drawnow

