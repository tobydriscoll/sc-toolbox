function [w,beta] = drawpoly(fig,cmd)
%DRAWPOLY Draw a polygon with the mouse.
%   P = DRAWPOLY allows the user to draw a polygon with the
%   mouse.  Use the mouse to position the crosshair and press the left
%   mouse button to create a vertex.  For use with other S-C Toolbox
%   functions, the vertices must be specified in a "positively oriented"
%   manner; i.e.  counterclockwise for interior polygons and clockwise
%   for exterior regions.  There are several GUI elements added to the
%   figure to help you snap vertices to a grid, get specfic angles, etc.
%   For the last vertex, use the middle or right mouse button, or double
%   click.  Upon return, P is a polygon object.
%
%   [W,BETA] = DRAWPOLY returns vectors of vertices and turning angles
%   instead of a polygon object.
%
%   See the user's guide for full details.
%	
%   See also POLYGON methods PLOT, EDIT, MODIFY.

%   Copyright 1998 by Toby Driscoll.
%   $Id: drawpoly.m 298 2009-09-15 14:36:37Z driscoll $

persistent youve_been_warned
if isempty(youve_been_warned)
  a = questdlg('DRAWPOLY is obsolete. Use POLYEDIT instead.','Obsolete function','Run DRAWPOLY','Run POLYEDIT','Run POLYEDIT');
  youve_been_warned = 1;
  if strcmp(a,'Run POLYEDIT')
    if nargin < 1
      w = polyedit;
    else
      w = polyedit(fig);
    end
    return
  end
end
  
% If called from the GUI, reroute to subfunction cback
if nargin > 1
  [w,beta] = cback(fig,cmd);
  return
end

if nargin < 1
  fig = gcf;
end

figure(fig);
turn_off_hold = ~ishold;

% Save current properties to be restored
figprops = {'unit';'pos';'resizefcn';'pointer';'windowbuttondown';...
    'windowbuttonup';'pointer'};
figstate = get(fig,figprops);
axprops = {'unit';'pos';'xgrid';'ygrid'};
axstate = get(gca,axprops);

% Get figure position and screen size
set(fig,'unit','point')
figpos = get(fig,'pos');
unit = get(0,'unit');
set(0,'unit','point')
maxpos = get(0,'screensize');
set(0,'unit',unit)

% Set up figure size
%figsize = [462 396];
%delta = figsize - figpos(3:4);
%figpos(1) = max(min(figpos(1),maxpos(3)-figsize(1)),0); % move left if needed
%figpos(2) = max(figpos(2)-delta(2),0); % move bottom down, not top up
%figpos(3:4) = figsize;
%set(fig,'resizefcn','')
%set(fig,'pos',figpos)
set(fig,'resizefcn','drawpoly(gcbo,''resize'');')

if ~ishold
  cla
  axis auto
end

% Reset axes limits?
if strcmp(get(gca,'xlimmode'),'auto') & strcmp(get(gca,'ylimmode'),'auto')
  axlim = [-4 4 -4 4];
else
  axlim = axis;
  axis(axlim)
end

view(2)
set(gca,'unit','point','box','on','xgrid','off','ygrid','off')
set(gca,'xlim',axlim(1:2),'ylim',axlim(3:4),'plotboxaspectratio',[1,1,1])
hold on

% Set up controls
dp_make_widgets(fig);
dp_resize_widgets(fig);

set(fig, 'pointer','crosshair');
preview = line(NaN,NaN,'linestyle','--','erasemode','xor',...
    'clipping','off','tag','preview');
vertices = line(NaN,NaN,'linesty','none','marker','.',...
    'markersize',12,'tag','vertices','erasemode','xor');
line(NaN,NaN,'visible','off','tag','DRP_PT');
  
% Preparation.
if ~strcmp(computer,'SUN4')
  % Kludge.  Draw preview line when button is pressed.
  set(fig,'windowbuttondownfcn','drawpoly(gcf,''move'');');
end
set(fig,'windowbuttonupfcn', 'drawpoly(gcf,''up'');');

set(preview,'user','first','xdata',[],'ydata',[]);
set(vertices,'xdata',[],'ydata',[],'zdata',[])
  
% Get vertices.
done = 0;
while ~done
  [x0,y0] = drawpoly(fig,'getpoint');
  x = get(vertices,'xdata');
  y = get(vertices,'ydata');
  atinf = get(vertices,'zdata');
  m = length(x);
  n = m - sum(atinf)/2;
  mode = get(preview,'user');
  axlim = axis;
  reflen = mean([diff(axlim(1:2)),diff(axlim(3:4))])/60;
  if m > 0 
    if norm([x0-x(1), y0-y(1)]) < reflen
      mode = 'done';
    end
  end
  x(m+1) = x0;  
  y(m+1) = y0;
  edges = findobj(fig,'tag','edges');
  
  switch mode
  case {'normal','first'}
    if n > 0
      edges(n) = plot(x(m:m+1),y(m:m+1),'-','erasemode','back');
      set(edges(n),'user',n)
    end
    if x0>=axlim(1) & x0<=axlim(2) & y0>=axlim(3) & y0<=axlim(4)
      % Finite vertex.
      atinf(m+1) = 0;
      set(preview,'user','normal')
    else
      atinf(m+1) = 1;
      set(preview,'user','infinite')
      set(fig,'pointer','cross')
    end
  case 'infinite'
    % Re-entry point for next edge.
    atinf(m+1) = 1;
    set(fig,'pointer','crosshair')
    set(preview,'user','return')
  case 'return'
    edges(n) = plot(x(m:m+1),y(m:m+1),'-','erasemode','back');
    set(edges(n),'user',n)
    atinf(m+1) = 0;
    set(preview,'user','normal')
  case 'done'
    x = x(1:m);
    y = y(1:m);
    done = 1;
  end  
  
  % Check buttonpress.
  if ~strcmp(get(fig,'selectiontype'),'normal') & (m > 0)
    done = 1;
  end
  
  set(vertices,'xdata',x,'ydata',y,'zdata',atinf);
  set(edges,'tag','edges')
  drawnow
end
  
m = length(x);
n = length(x) - sum(atinf)/2;
if n > 0
  edges(n) = plot(x([m,1]),y([m,1]),'-');
end
drawnow

w = x(:)+i*y(:);
beta = scangle(w);
s = sum(atinf);
if s > 0
  b = beta(logical(atinf));
  b(b>0) = b(b>0) - 2;
  b = sum(reshape(b,2,s/2))';
  idx = find(atinf);
  w(idx(1:2:s)) = Inf*ones(s/2,1);
  beta(idx(1:2:s)) = b;
  w(idx(2:2:s)) = [];
  beta(idx(2:2:s)) = [];
end

% Clean up the mess.
set(fig,figprops,figstate)
set(gca,axprops,axstate)
set(vertices,'erasemode','normal')
set(edges,'erasemode','normal')
frame(1) = findobj(fig,'tag','dp_settings');
frame(2) = findobj(fig,'tag','dp_actions');
delete(get(frame(1),'userdata'))
delete(get(frame(2),'userdata'))
delete(frame)
delete(preview)
delete(vertices)
delete(edges)
  
% Construct a polygon object
if sum(beta) < 0
  poly = polygon(w,1+beta);
else
  poly = polygon(flipud(w),1-flipud(beta));
end
if nargout < 2
  w = poly;
end

axis auto
if turn_off_hold
  hold off
end
plot(poly)
set(gca,'xtickmode','auto','ytickmode','auto')
set(gca,'xticklabelmode','auto','yticklabelmode','auto')


function [out1,out2] = cback(fig,cmd)

% Invocations as a callback
if nargout > 0
  out1 = [];
  if nargout > 1
    out2 = [];
  end
end

switch cmd

case 'getpoint'
  DRP_PT = findobj(fig,'tag','DRP_PT');
  set(fig,'windowbuttonmotionfcn','drawpoly(gcf,''move'');');
  drawnow
  waitfor(DRP_PT,'user')
  data = get(DRP_PT,'user');
  out1 = data(1);
  out2 = data(2);
  
case 'move'
  ptrpos = get(gca,'currentpoint');
  P = ptrpos(1,1:2);
  axlim = axis;
  preview = findobj(fig,'tag','preview');
  vertices = findobj(fig,'tag','vertices');
  x = get(vertices,'xdata');
  y = get(vertices,'ydata');
  pts = [x(:),y(:)];
  m = size(pts,1);
  mode = get(preview,'user');
  grid = str2num(get(findobj(fig,'tag','grid_space'),'string'));
  grid = grid * get(findobj(fig,'tag','snap'),'value');
  qang = 1/str2num(get(findobj(fig,'tag','ang_quant'),'string'));
  qang = qang * get(findobj(fig,'tag','qang'),'value');
  qlen = str2num(get(findobj(fig,'tag','len_quant'),'string'));
  qlen = qlen * get(findobj(fig,'tag','qlen'),'value');
  
  % Modify point to meet mode constraints.
  if strcmp(mode,'normal')
    % No constraints.
    mode = 0;
  elseif strcmp(mode,'first')
    % Point may not be outside axes box.
    P(1) = min(max(P(1),axlim(1)),axlim(2));
    P(2) = min(max(P(2),axlim(3)),axlim(4));
    mode = 0;
  elseif strcmp(mode,'infinite')
    % Point may not be inside axes box.
    if all(P>axlim([1,3])) & all(P<axlim([2,4]))
      [junk,j] = min(abs([P(1)-axlim(1:2);P(2)-axlim(3:4)]'));
      P = axlim(j);
    end
    mode = 1;
  elseif strcmp(mode,'return')
    % Point may not be outside axes box.
    P(1) = min(max(P(1),axlim(1)),axlim(2));
    P(2) = min(max(P(2),axlim(3)),axlim(4));
    ang = scangle([pts(m-2:m,1);P(1)]+i*[pts(m-2:m,2);P(2)]);
    ang = ang(2:3);
    ang(ang>0) = ang(ang>0) - 2;
    ang = sum(ang);
    mode = 2;
  end

  % Modify point to meet angle, length, or grid constraints.
  
  if any(mode==[0,2]) & qang & (m > 0)	% quantized angle
    % Find arg of new side which meets quantization requirements. 
    if m==1
      % No reference side, so refer to positive real axis
      ang = angle(P(1)-pts(m,1)+i*(P(2)-pts(m,2)))/pi;
      theta = qang*round(ang/qang)*pi;
    elseif mode==0
      ang = scangle([pts(m-1:m,1);P(1)]+i*[pts(m-1:m,2);P(2)]);
      ang = qang*round(ang(2)/qang);
      theta = atan2(pts(m,2)-pts(m-1,2),pts(m,1)-pts(m-1,1))-pi*ang;
    elseif mode==2
      % ang was computed above
      ang = qang*round(ang/qang);
      theta = atan2(pts(m-1,2)-pts(m-2,2),pts(m-1,1)-pts(m-2,1))-pi*ang;
    end
    % Project P to correct angle.
    A = pts(m,:);
    BA = [cos(theta),sin(theta)];
    P = A + ((BA)*(P-A)')*(BA);
    grid = 0;
  end
  if (mode==0) & qlen & (m > 0)		% quantized length
    A = pts(m,:);
    len = norm(P-A);
    fixlen = qlen*(round(len/qlen));
    P = A + fixlen/len*(P-A);
    grid = 0;
  end
  if any(mode==[0,1,2]) & grid	        % snap to grid
    grid = [grid,grid];
    minxy = axlim([1,3]);
    P = minxy + grid.*(round((P-minxy)./grid));
  end
  
  % If returning from infinity, don't allow angle > -1.
  if mode==2
    if ang > -1+2*eps
      % Would be illegal.  Project to make ang=-1.
      A = pts(m,:);
      B = pts(m-2,:) + A - pts(m-1,:);
      P = A + ((B-A)*(P-A)')/((B-A)*(B-A)')*(B-A);
      ang = -1;
      qang = 0;				% override other restrictions
      grid = 0;
    elseif ang < -3
      % It's illegal.  Is it even possible?
      P = [NaN,NaN];
      ang = NaN;
    end
  end

  % Update.
  if m > 0
    if ~(mode==1)			% preview line
      set(preview, 'xdata',[pts(m,1),P(1)], 'ydata',[pts(m,2),P(2)]);
    end
  end
  set(vertices,'userdata',P);
  drawnow

case 'up'	
  preview = findobj(fig,'tag','preview');
  vertices = findobj(fig,'tag','vertices');
  P = get(vertices,'userdata');
  if ~isnan(P(1))			% valid point
    set(gcf,'windowbuttonmotionfcn','');
    set(preview,'xdata',[P(1),NaN], 'ydata',[P(2),NaN])
    set(vertices,'userdata',[NaN,NaN])
    set(findobj(fig,'tag','DRP_PT'),'user',P);
  end

case 'snap'  
  snap = findobj(fig,'tag','snap');
  qang = findobj(fig,'tag','qang');
  qlen = findobj(fig,'tag','qlen');
  if get(snap,'value')
    set(qang,'value',0)
    set(qlen,'value',0)
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

case 'qang'  
  snap = findobj(fig,'tag','snap');
  qang = findobj(fig,'tag','qang');
  if get(qang,'value')
    set(snap,'value',0)
    drawpoly(fig,'snap');
  end

case 'qlen'  
  snap = findobj(fig,'tag','snap');
  qlen = findobj(fig,'tag','qlen');
  if get(qlen,'value')
    set(snap,'value',0)
    drawpoly(fig,'snap');
  end
  
case 'xlim'
  xlim = str2num(get(findobj(fig,'tag','xlimits'),'string'));
  set(gca,'xlim',xlim)
  drawpoly(fig,'snap');
  
case 'ylim'
  ylim = str2num(get(findobj(fig,'tag','ylimits'),'string'));
  set(gca,'ylim',ylim)
  drawpoly(fig,'snap');
  
case 'clear'
  preview = findobj(fig,'tag','preview');
  vertices = findobj(fig,'tag','vertices');
  set(preview,'user','first','xdata',[NaN NaN],'ydata',[NaN NaN]);
  set(vertices,'xdata',[],'ydata',[],'zdata',[])
  delete(findobj(fig,'tag','edges'))
  set(fig,'color',get(fig,'color'))
  
case 'undo'
  preview = findobj(fig,'tag','preview');
  vertices = findobj(fig,'tag','vertices');
  x = get(vertices,'xdata');
  y = get(vertices,'ydata');
  atinf = get(vertices,'zdata');
  mode = get(preview,'user');
  P = get(vertices,'user');
  m = length(x);
  n = m - floor(sum(atinf)/2);
  mask = ones(m,1);
  kill_edge = 0;
  if strcmp(mode,'normal') & (m > 0)
    mask(m) = 0;
    if m > 1
      if atinf(m-1)
	set(preview,'user','return')
      end
    else
      set(preview,'user','first')
    end
    kill_edge = 1;
  elseif strcmp(mode,'infinite')
    mask(m) = 0;
    set(preview,'user','normal')
    set(fig,'pointer','crosshair')
    kill_edge = 1;
  elseif strcmp(mode,'return')
    mask(m) = 0;
    set(preview,'user','infinite')
    set(preview,'xdata',[NaN,NaN],'ydata',[NaN,NaN])
    set(fig,'pointer','cross')
  end
  mask = logical(mask);
  set(vertices,'erasemode','norm')
  set(vertices,'xdata',x(mask),'ydata',y(mask),'zdata',atinf(mask));
  set(vertices,'erasemode','back')
  edges = findobj(fig,'tag','edges');
  if kill_edge & (n > 1)
    delete(findobj(edges,'user',n-1))
  end
  if ~any(mask)				% no vertices left
    % Erase preview line 
    set(findobj(fig,'tag','preview'),'xdata',NaN,'ydata',NaN,'user','first')
  end
  set(gcf,'color',get(gcf,'color'))
  drawpoly(fig,'move');
  
case 'done'
  preview = findobj(fig,'tag','preview');
  vertices = findobj(fig,'tag','vertices');
  set(preview,'user','done')
  set(fig,'windowbuttonmotion','')
  P = get(vertices,'userdata');
  set(preview,'xdata',[P(1),NaN], 'ydata',[P(2),NaN])
  set(findobj(fig,'tag','DRP_PT'),'user',[NaN,NaN]);
  
case 'cancel'
  preview = findobj(fig,'tag','preview');
  vertices = findobj(fig,'tag','vertices');
  set(preview,'user','done')
  set(fig,'windowbuttonmotion','')
  set(vertices,'xdata',[],'ydata',[],'zdata',[])
  set(findobj(fig,'tag','DRP_PT'),'user',[NaN,NaN]);
  
case 'resize'
  dp_resize_widgets(fig);
  
end
  

function dp_make_widgets(fig)

bgcolor = .7*[1 1 1];
set(fig,'defaultuicontrolbackground',bgcolor,'defaultuicontrolunits','point')

% We aren't going to set positions, because that is done in dp_resize_widgets

% Create uicontrol frames
frame(1) = uicontrol('style','frame','tag','dp_settings');
frame(2) = uicontrol('style','frame','tag','dp_actions');
  
widget(1) = uicontrol('style','radio',...
    'string','Snap to grid',...
    'call','drawpoly(gcf,''snap'');',...
    'horizontal','l',...
    'tag','snap');
widget(2) = uicontrol('style','radio',...
    'string','Quantize angle',...
    'call','drawpoly(gcf,''qang'');',...
    'hor','l',...
    'tag','qang');
widget(3) = uicontrol('style','radio',...
    'string','Quantize length',...
    'call','drawpoly(gcf,''qlen'');',...
    'hor','l',...
    'tag','qlen');

widget(4) = uicontrol('style','text',...
    'string','Grid spacing:',...
    'hor','r');
widget(5) = uicontrol('style','text',...
    'string','Angle quantum:',...
    'hor','r');
widget(6) = uicontrol('style','text',...
    'string','Length quantum:',...
    'hor','r');

gridspace = diff(get(gca,'xlim'))/16;
widget(7) = uicontrol('style','edit',...
    'string',sprintf('%.3g',gridspace),...
    'background',.8*[1 1 1],...
    'call','drawpoly(gcf,''snap'');',...
    'tag','grid_space');
widget(8) = uicontrol('style','text',...
    'string','p / ',...
    'fontname','symbol',...
    'fontsize',11,...
    'hor','r');
widget(9) = uicontrol('style','edit',...
    'string','6',...
    'background',.8*[1 1 1],...
    'tag','ang_quant');
widget(10) = uicontrol('style','edit',...
    'string',sprintf('%.3g',2*gridspace),...
    'background',.8*[1 1 1],...
    'tag','len_quant');

axlim = axis;
widget(11) = uicontrol('style','text',...
    'string','x limits:',...
    'hor','r');
widget(12) = uicontrol('style','text',...
    'string','y limits:',...
    'hor','r');
widget(13) = uicontrol('style','edit',...
    'string',sprintf('[%.3g %.3g]',axlim(1:2)),...
    'background',.8*[1 1 1],...
    'tag','xlimits',...
    'call','drawpoly(gcf,''xlim'');');
widget(14) = uicontrol('style','edit',...
    'string',sprintf('[%.3g %.3g]',axlim(3:4)),...
    'background',.8*[1 1 1],...
    'tag','ylimits',...
    'call','drawpoly(gcf,''ylim'');');

widget(15) = uicontrol('style','push',...
    'string','Clear',...
    'call','drawpoly(gcf,''clear'');');
widget(16) = uicontrol('style','push',...
    'string','Undo',...
    'call','drawpoly(gcf,''undo'');');
widget(17) = uicontrol('style','push',...
    'string','Done',...
    'call','drawpoly(gcf,''done'');');
widget(18) = uicontrol('style','push',...
    'string','Cancel',...
    'call','drawpoly(gcf,''cancel'');');

set(frame(1),'userdata',widget(1:14))
set(frame(2),'userdata',widget(15:18))


function dp_resize_widgets(fig)

figpos = get(fig,'pos');

% Frames
f(1) = findobj(fig,'tag','dp_settings');
f(2) = findobj(fig,'tag','dp_actions');

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
[2 49 100 20];...
[2 26 100 20];...
[2  3 100 20];...
[120 52 80 14];...
[120 28 80 14];...
[120 4 80 14];...
[204 50 48 18];...
[204 27 20 14];...
[224 27 28 18];...
[204 4 48 18];...
[268 52 48 14];...
[268 31 48 14];...
[320 51 54 18];...
[320 30 54 18]
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
