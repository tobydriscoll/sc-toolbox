function out = scgui(varargin)
%SCGUI  Create graphical user interface for the SC Toolbox.
%   SCGUI creates the graphical user interface (GUI) for the
%   Schwarz-Christoffel Toolbox in a new figure window.
%
%   For complete details on the interface, see the user's guide.

%   Copyright 1998-2002 by Toby Driscoll.
%   $Id: scgui.m 298 2009-09-15 14:36:37Z driscoll $

% If called from the GUI, reroute to switchyard function cback
%%if nargin > 1
%%  out = cback(guidata(gcbo),varargin{:});
%%  return
%%end

% Create figure
fig = figure('name','Schwarz-Christoffel Mapping','numbertitle','off',...
    'tag','scgui','menubar','none','vis','off');
set(fig,'defaultuicontrolinterrupt','on')

% Determine/set figure size
minsize = [95 36];
unit = get(0,'units');
set(0,'units','char')
ss = get(0,'screensize');
set(0,'units',unit)
optsize = 0.65*ss(3:4);
truesize = max(optsize,minsize);
set(fig,'unit','char','pos',[0 ss(4)-truesize(2)-8 truesize])
% can't figure out why this doesn't work right.
%% movegui(fig,'northeast')

% Custom resize function
% I have to disable this because at least in my window manager (KDE), the
% resizing is soooooooo slow as to be painful. But you can try it on
% your machine.
%%set(fig,'resizefcn',@resize,'deletefcn',@selfdelete)
set(fig,'resize','off')

% Create axes
ax(1) = axes('units','char','tag','PhysicalAxes');
ax(2) = axes('units','char','tag','CanonicalAxes');
set(ax,'next','add','box','on','plotboxaspectratio',[1 1 1])

% Create widgets and set sizes appropriately (subfunctions)
make_widgets(fig);
resize_widgets(fig)
drawnow

set(fig,'closerequestfcn',@scquit)
set(fig,'vis','on')


function draw(obj,varargin)
data = guidata(obj);

p = polyedit;

if ~isempty(p)
  deletepoly(obj,'confirmed')
  deletemap(obj)
  addpoly(obj,p)
end


function domodify(obj,varargin)
data = guidata(obj);
p = data.polygon;
if isempty(p)
  errordlg('You must first draw or import a polygon,','SC Error')
  return
end		
n = length(p);

p1 = polyedit(p);

% Was anything really changed?
changed = 0;
w = vertex(p);
w1 = vertex(p1);
if length(w1) ~= length(w)
  changed = 1;
elseif norm(w(~isinf(w))-w1(~isinf(w1))) > 1000*eps
  changed = 1;
end

if changed
  deletepoly(obj,'confirmed')
  addpoly(obj,p1)
  % Ask about continuation.
  set(findobj(data.SCfig,'tag','MeshLines'),'erase','norm','vis','off')
  if data.iscurrent & length(w1)==length(w)
    a = questdlg('Do you want to use continuation from the current parameter solution?','Map continuation','Yes','No','No');
  else
    a = 'No';
  end
  if strcmp(a,'No')
    deletemap(obj)
  else
    deletemap(obj,'continue')
  end
end



function edit(obj,varargin)
data = guidata(obj);
w = vertex(data.polygon);
if isempty(w), return, end
if any(isinf(w))
  errordlg('Cannot edit vertices of unbounded polygons.','SC Error')
  return
end

n = length(w);

% Call edit window
[flag,w] = scgedit('Edit vertices','Vertices:',w);
% If edit was accepted, try to convert it to a polygon.
if flag > 0
  try 
    p = polygon(w);
  catch
    errordlg('Invalid polygon specified.','SC Error')
    return
  end
else
  return
end

% We now assume the polygon was changed.
deletepoly(obj,'confirmed')
if length(w)~=n
  % No way to be a continuation
  deletemap(obj)
else
  deletemap(obj,'continue')
end

addpoly(obj,polygon(w))


function solve(obj,varargin)
data = guidata(obj);
p = data.polygon;
if isempty(p)
  errordlg('You must first draw or import a polygon.','SC Error')
  return
end

mapnum = get(data.DomainPopup,'value');
maptype = data.mapclass{mapnum};
oldmap = data.map;
trace = 1;
%trace = get(data.TraceBox,'value');
tol = str2num(get(data.TolEdit,'string'));

if trace == 1
  trace = 'on';
else
  trace = 'off';
end

opt = sctool.scmapopt('Trace',trace,'Tol',tol);

if ~isempty(oldmap)
  % Continuation
  mapcmd = sprintf('map = %s(oldmap,p,opt);',maptype);
else
  mapcmd = sprintf('map = %s(p,opt);',maptype);
end

try
  eval(mapcmd)
catch
  errordlg('Solution process failed!','SC Error')
  return
end
    
deletemap(obj)
% Store new map
if ~isempty(map)
  addmap(obj,map)
  % Get the polygon back out of the map, because indexing and trivial
  % vertices may change.
  p = polygon(map);
  deletepoly(obj,'confirmed')
  addpoly(obj,p)
end

plotcanonical(obj)
setview(obj)


function plotcanonical(obj,varargin)
data = guidata(obj);
fig = data.SCfig;
if ~data.iscurrent, return, end

% Draw the canonical domain
axes(data.CanonicalAxes)
delete(findobj(gca,'tag','PolygonPlot'))
z = parameters(data.map);
z = z.prevertex;

switch class(data.map)
  case { 'diskmap','extermap','crdiskmap' }
    h = plot(exp(1i*linspace(0,2*pi,100)),'linewid',1);
    axis normal
    axis equal
    axis square
    axis([-1.1 1.1 -1.1 1.1])
  case 'hplmap'
    xmin = min(z(~isinf(z)));
    xmax = max(z(~isinf(z)));
    xlim = [1.1*xmin-.1*xmax,1.1*xmax-.1*xmin];
    h = plot(xlim,[0,0],'linewid',1);
    axis normal
    axis equal
    axis square
    axis([xlim,-.1,1])
  case 'stripmap'
    minx = min(real(z(~isinf(z))));
    maxx = max(real(z(~isinf(z))));
    d = (maxx-minx)/2;
    m = (maxx+minx)/2;
    h(1) = plot([m-1.2*d,m+1.2*d],[0,0],'linewid',1);
    h(2) = plot([m-1.2*d,m+1.2*d],[1,1],'linewid',1);
    axis normal
    axis equal
    axis([m-1.2*d,m+1.2*d,-.1,1.1])
  case { 'rectmap','crrectmap' }
    axis auto
    h = plot(polygon(z));
    axis equal
end 	

guidata(fig,data)
    
set(h,'tag','PolygonPlot')
h = plot(real(z),imag(z),'ko','markersize',5,'markerfacecolor','k');
set(h,'tag','PolygonPlot')


function recenter(obj,varargin)
data = guidata(obj);
fig = data.SCfig;
% Only applies to current disk maps
isdisk = any(strcmp(class(data.map),{'diskmap','crdiskmap'}));
if ~data.iscurrent
  errordlg('You must first solve for the map parameters.','SC Error')
  return
elseif ~isdisk
  errordlg('The conformal center is defined only for disk maps.','SC Error')
  return
end

ax(1) = data.PhysicalAxes;
ax(2) = data.CanonicalAxes;
% Show physical domain & turn off point mapping
set(findobj(ax(1)),'vis','on')
set(findobj(ax(2)),'vis','off')
set(ax,'buttondown','')

% Get point and do calculation
axes(ax(1))
title('Click at conformal center')
[xc,yc] = ginput(1);
title('')
data.map = center(data.map,xc+i*yc);
guidata(obj,data)

deletemap(obj)
addmap(obj,data.map)

plotcanonical(obj)
setview(obj)

  
function scdisplay(obj,varargin)
data = guidata(obj);
if ~data.iscurrent
  errordlg('You must first solve for the map parameters.','SC Error')
else
  str = char(data.map);
  pos = [0 0 100 size(str,1)];
  f = figure('unit','char','vis','off','integerhandle','off',...
      'name','Map parameters','numbertitle','off','pos',pos);
  u = uicontrol('style','text','unit','char',...
      'pos',pos,'fontname','FixedWidth','horiz','left',...
      'background','white','string',str);
  set(f,'vis','on')
end


function meshplot(obj,varargin)
data = guidata(obj);
if ~data.iscurrent
  errordlg('You must first solve for the map parameters.','SC Error')
  return
end

ax(1) = data.PhysicalAxes;
ax(2) = data.CanonicalAxes;

% Parse mesh line widgets
xrhan = data.XREdit;
ythan = data.YTEdit;
xrstring = get(xrhan,'string');
ytstring = get(ythan,'string');
if strcmp(xrstring,'defaults') | strcmp(xrstring,'default')
  xr = 10;
else
  xr = str2num(xrstring);
end
if strcmp(ytstring,'defaults') | strcmp(ytstring,'default')
  yt = 10;
else
  yt = str2num(ytstring);
end

% Bring up physical domain and clear old mesh lines
%%    set(findobj(ax(1)),'vis','on')
%%    set(findobj(ax(2)),'vis','off')
ml = findobj(ax,'tag','MeshLines');
%%set(ml,'erase','normal')
delete(ml)
delete(findobj(ax(1),'tag','PolygonPlot'))
drawnow

% Exterior map: Allow axis autoscale
if strcmp(class(data.map),'extermap')
  set(ax(1),'xlimmode','auto','ylimmode','auto')
end

% Do the plot
[h,xr,yt] = plot(data.map,xr,yt);
%%set(h,'erasemode','none','tag','MeshLines')
set(h,'tag','MeshLines')

% Plotting routines plot their own copy of the polygon,
% without returning handles. This locates and (re)tags these lines.
hp = findobj(ax(1),'tag','PolygonSide');
set(hp,'tag','PolygonPlot','uicontext',data.PolygonContextMenu)

% Update mesh line widgets
if isempty(xr)
  xrstring = '';
else
  xrstring = ['[',sprintf('%.3g ',xr),']'];
end
isdisk = strcmp(class(data.map),'diskmap');
isdisk = isdisk | strcmp(class(data.map),'extermap');    
if isempty(yt)
  ytstring = '';
elseif isdisk
  ytstring = ['pi*[',sprintf('%.3g ',yt/pi),']'];
else
  ytstring = ['[',sprintf('%.3g ',yt),']'];
end
set(xrhan,'string',xrstring)
set(ythan,'string',ytstring)

% May want to fix up the canonical domain
axes(ax(2))
z = prevertex(data.map);
switch class(data.map)
  case {'diskmap','extermap','crdiskmap'}
  case 'hplmap'
    zf = z(~isinf(z));
    xlim = [min(zf);max(zf)];
    xlim = [1.1 -.1; -.1 1.1] * xlim;
    ylim = [0;max(zf)];
    ylim = [0 0; -.1 1.1]*ylim;
    axis([xlim',-0.05*ylim(2),ylim(2)])
  case 'stripmap'
    %%	zf = z(~isinf(z));
    %%	xlim = [min(zf);max(zf)];
    %%	xlim = [1.1 -.1; -.1 1.1] * xlim;
    %%	axis([xlim',-.1,1.1])
  case 'rectmap'
  case 'crrectmap'
end
% Want to do this so that the buttondown of the mesh lines is set up.
setview(obj)


function images = mappts(pts)
domain = get(gca,'tag');
data = guidata(gcbo);
source = pts;				% (given)
images = [];
if data.iscurrent
  title('Mapping...')
  drawnow
  
  if ~isempty(source)
    % Forward/inverse map, depending on plane in view
    if strcmp(domain,'PhysicalAxes')
      images = evalinv(data.map,source);
    elseif strcmp(domain,'CanonicalAxes')
      images = eval(data.map,source);
    end
  else
    images = [];
  end

  title('')
  drawnow
end


function addpt(varargin)
% Take no action unless it's a standard left-click. But DON'T just return,
% because that short-circuits context menu clicks.
if strcmp(get(gcf,'selectiontype'),'normal')
  data = guidata(gcbo);
  pt = get(gca,'currentpoint');
  % Map it
  img = mappts(pt(1,1)+i*pt(1,2));
  ax(1) = data.PhysicalAxes;
  ax(2) = data.CanonicalAxes;
  domain = find(gca==ax);

  % Ensure the existence of the line objects.
  if isempty(data.PhysicalPoints) | ~ishandle(data.PhysicalPoints)
    data.PhysicalPoints = make_pointsobject(data,'Physical');
  end
  if isempty(data.CanonicalPoints) | ~ishandle(data.CanonicalPoints)
    data.CanonicalPoints = make_pointsobject(data,'Canonical');
  end
  h = [data.PhysicalPoints data.CanonicalPoints];
  
  % Update objects to show points
  axes(ax(domain))
  %%set(h,'erasemode','none')
  xs = [get(h(domain),'xdata') pt(1,1)];
  ys = [get(h(domain),'ydata') pt(1,2)];
  set(h(domain),'xdata',xs,'ydata',ys)
  xi = [get(h(3-domain),'xdata') real(img)];
  yi = [get(h(3-domain),'ydata') imag(img)];
  set(h(3-domain),'xdata',xi,'ydata',yi)
  if domain==1
    data.phypoints = xs + 1i*ys;
    data.canpoints = xi + 1i*yi;
  else
    data.phypoints = xi + 1i*yi;
    data.canpoints = xs + 1i*ys;
  end
  %%set(h,'erasemode','normal')
  
  guidata(gcbf,data)
end


function editpts(varargin)
domain = get(gca,'tag');
if isempty(domain), return, end
data = guidata(gcbo);
pphan = data.PhysicalPoints;
cphan = data.CanonicalPoints;

wp = get(pphan,'xdata')+i*get(pphan,'ydata');
zp = get(cphan,'xdata')+i*get(cphan,'ydata');
wp = wp(~isnan(wp));
zp = zp(~isnan(zp));

% Interpretation of box depends on current view domain
if strcmp(domain,'PhysicalAxes')
  [flag,wp] = scgedit('Edit points','Physical',wp);
  if flag > 0
    zp = mappts(wp);
  end
elseif strcmp(domain,'CanonicalAxes')
  [flag,zp] = scgedit('Edit points','Canonical',zp);
  if flag > 0
    wp = mappts(zp);
  end
end

if flag > 0
  % Update point objects
  set(pphan,'xdata',real(wp(:))','ydata',imag(wp(:))')
  set(cphan,'xdata',real(zp(:))','ydata',imag(zp(:))')
  data.phypoints = wp(:);
  data.canpoints = zp(:);
  guidata(gcbf,data)
end


function importexport(obj,varargin)
data = guidata(obj);
fig = data.SCfig;
% Find import/export window, or create if necessary
fig2 = findobj(0,'tag','sc_importexport');
if isempty(fig2)
  fig2 = make_importexportfig;
else
  figure(fig2)
end					% (create window)
data = guidata(fig2);
data.parentfig = fig;
guidata(fig2,data)

% Reveal window and wait for action
set(fig2,'vis','on')
uiwait(fig2)
set(fig2,'vis','off')			% turn it off


function scimport(obj,varargin)
% Data from the I/E dialog.
IEdata = guidata(obj);
% Data from the SC figure
SCfig = IEdata.parentfig;
SCdata = guidata(SCfig);

name = get(IEdata.PolygonEdit,'string');
if ~isempty(name)
  try 
    p = evalin('base',name);
    deletepoly(SCfig,'confirmed')
    deletemap(SCfig)
    addpoly(SCfig,p)
  catch
    errordlg(lasterr,'SC Error')
  end
end

name = get(IEdata.MapEdit,'string');
if ~isempty(name)
  try 
    map = evalin('base',name);
    % This overrides an imported polygon.
    p = polygon(map);
    deletepoly(SCfig,'confirmed')
    deletemap(SCfig)
    addpoly(SCfig,p)
    addmap(SCfig,map)
  catch
    errordlg(lasterr,'SC Error')
    uiresume(gcbf)
  end
end

plotcanonical(SCfig)
setdomain(SCfig)
setview(SCfig)
uiresume(gcbf)


function scexport(varargin)
% Data from the I/E dialog.
IEdata = guidata(gcbf);
% Data from the SC figure
SCfig = IEdata.parentfig;
SCdata = guidata(SCfig);

name = get(IEdata.PolygonEdit,'string');
if ~isempty(name) & ~isempty(SCdata.polygon)
  assignin('base',name,SCdata.polygon);
end
  
name = get(IEdata.MapEdit,'string');
if ~isempty(name) & SCdata.iscurrent
  assignin('base',name,SCdata.map);
end

uiresume(gcbf)


function setview(obj,varargin)
data = guidata(obj);
ax(1) = data.PhysicalAxes;
ax(2) = data.CanonicalAxes;
viewdom = get(data.ViewPopup,'value');

% Set axes positions
window = data.axeswindow;
if viewdom==1 | viewdom==2
  set(ax,'pos',window)
else
  % They must occupy the window together
  dx = .46*window(3);
  dy = window(4);
  set(ax(1),'pos',[window(1) window(2) dx dy])
  set(ax(2),'pos',[window(1)+window(3)-dx window(2) dx dy])
end

% Set visibility, map-on-click functions
vis = [viewdom~=2 viewdom~=1];
set(findobj(ax(vis)),'visible','on')
set(findobj(ax(~vis)),'visible','off')
% All children get the buttondown set to map points. This is because lines
% can "hide" the axes in a zone around the line (so no mapping close to the
% boundary). Even though the 'hittest' property seems to fix this, 
% if that is turned off then we don't get nice context menus. The addpt
% function will sort out left clicks from right.
if ~isempty(data) & data.iscurrent
  set(findobj(ax(vis)),'buttondown',@addpt)
else
  set(findobj(ax(vis)),'buttondown','')
end
set(findobj(ax(~vis)),'buttondown','')

if sum(vis)==1
  axes(ax(vis))
end


function setdomain(obj,varargin)
data = guidata(obj);
fig = data.SCfig;
mapnum = get(data.DomainPopup,'value');
%type = {'disk','half-plane','strip','rectangle','exterior'};
curmapnum = strmatch(class(data.map),data.mapclass);
  
if ~isempty(curmapnum) & mapnum~=curmapnum
  % Map type has been changed, so kill old map
  deletemap(obj)
end
  
% Mesh line options should be r/theta or x/y
set(data.XREdit,'string','defaults');
set(data.YTEdit,'string','defaults');
xrtext = data.XRText;
yttext = data.YTText;
if mapnum==1 | mapnum==5
  set(xrtext,'string','r =')
  set(yttext,'string','theta =')
else
  set(xrtext,'string','x =')
  set(yttext,'string','y =')
end
  
%%% Clear all current mesh lines 
%%delete(findobj(fig,'tag','MeshLines'))


function scquit(varargin)
if strcmp(questdlg('Really quit?','Quit SC mapping'),'Yes')
  selfdelete
end


function selfdelete(varargin)
delete([findobj(0,'tag','sc_importexport') gcbf])


function resize(varargin)
resize_widgets(gcbf)


function savedata(varargin)
obj = gcbo;
data = guidata(obj);
p = data.polygon;  map = data.map;
if isempty(p) && isempty(map)
  return
end
if ~isfield(data,'savename') || isempty(data.savename)
  filter = '*.mat';
else
  filter = data.savename;
end
[filen,pathn] = uiputfile(filter,'Save SC data');
if isstr(filen)
  data.savename = [pathn filen];
  guidata(obj,data);
  save([pathn filen],'p','map');
end


function loaddata(varargin)
obj = gcbo;
data = guidata(obj);
if ~isfield(data,'savename') || isempty(data.savename)
  filter = '*.mat';
else
  filter = data.savename;
end
[filen,pathn] = uigetfile(filter,'Load SC data');
if isstr(filen)
  S = load([pathn filen]);
  if ~isfield(S,'p') || ~isfield(S,'map')
    error('File %s does not contain needed GUI data',[pathn filen])
  else
    data.savename = [pathn filen];
    guidata(obj,data);
    deletepoly(obj,'confirmed')
    deletemap(obj)
    addpoly(obj,S.p)
    if ~isempty(S.map), addmap(obj,S.map), end
  end
end


function printablecopy(varargin)
fig = gcbf;
data = guidata(fig);
ax(1) = data.PhysicalAxes;
ax(2) = data.CanonicalAxes;
viewdom = get(data.ViewPopup,'value');

if viewdom==1 || viewdom==2  % only one view
  idx = viewdom;
else
  idx = [1 2];
end

newfig = figure;
for l = 1:length(idx)
  subplot(1,length(idx),l);
  tmp = gca;
  axun = get(tmp,'unit');  axpos = get(tmp,'pos');
  newax = copyobj( ax(idx(l)), newfig );
  delete(tmp)
  set(newax,'unit',axun,'pos',axpos);
  % Wipe out GUI-specific markers.
  set(newax,'buttondown','','tag','')
  set(allchild(newax),'uicontextmenu',[],'buttondown','','tag','')
end
  

function addpoly(obj,p,varargin)
% Add a polygon. Make sure that it is displayed and tagged properly, and
% enable appropriate further actions.
data = guidata(obj);
data.polygon = p;
guidata(obj,data)

if isempty(p), return, end

if nargin<3 | ~all(ishandle(varargin{1}))
  axes(data.PhysicalAxes)
  axis auto
  h = plot(p);
else
  h = varargin{1};
end
set(h, 'tag','PolygonPlot', 'uicontext',data.PolygonContextMenu);

set(data.PolygonEnabled,'enable','on')


function deletepoly(obj,varargin)
% Remove the polygon from display and the data structure.
if nargin < 2 | ~isequal( varargin{1},'confirmed')
  answer = questdlg('Really delete the current polygon?','Confirm delete',...
      'Yes','No','No');
  if strcmp(answer,'No'), return, end
end

data = guidata(obj);
data.polygon = [];
guidata(obj,data)

delete(findobj(data.SCfig,'tag','PolygonPlot'))

% Without a polygon, some controls are disabled.
set(data.PolygonEnabled,'enable','off')


function menudeletepoly(obj,varargin)
deletepoly(obj)
deletemap(obj)


function addmap(obj,map)
% Add a map. Make sure that supporting objects exist and 
% enable appropriate further actions.
data = guidata(obj);
data.map = map;
data.iscurrent = 1;
guidata(obj,data)

% Correct the class popup to reflect reality.
set(data.DomainPopup,'value',strmatch(class(map),data.mapclass))
plotcanonical(obj)

% Enable relevant controls.
set(data.MapEnabled,'enable','on')


function deletemap(obj,varargin)
% The current map has been deleted (or never existed). Take care of the data
% structures and widgets.
data = guidata(obj);
% We may not erase the memory of the map if we want to do continuation.
if nargin < 2 | ~isequal(varargin{1},'continue')
  data.map = [];
end
data.iscurrent = 0;
guidata(obj,data)

set(data.MapEnabled,'enable','off')

% No mesh lines or canonical domain.
ml = findobj(data.SCfig,'tag','MeshLines');
%%set(ml,'erasemode','normal')
delete(ml)
delete(findobj(data.CanonicalAxes,'tag','PolygonPlot'))

% There can be no points to display...
deletepoints(obj,'confirmed')
% ...nor can more be mapped.
set(findobj([data.PhysicalAxes,data.CanonicalAxes]),'buttondown','')


function deletepoints(obj,varargin)
% Remove plotted points from display and data structure.
if nargin < 2 | ~isequal( varargin{1},'confirmed')
  answer = questdlg('Really delete the plotted points?','Confirm delete',...
      'Yes','No','No');
  if strcmp(answer,'No'), return, end
end

data = guidata(obj);
data.phypoints = [];
data.canpoints = [];
h = [data.PhysicalPoints,data.CanonicalPoints];
set(h(ishandle(h)),'xdata',[],'ydata',[])
guidata(obj,data)



function hp = make_pointsobject(data,type)
hp = line('parent',getfield(data,[type 'Axes']),...
    'linestyle','none','marker','o','markersize',5,...
    'color',[1 0 0],'markerfacecolor',[1 0 0],...
    'tag',[type 'Points'],'uicontextmenu',data.PointContextMenu,...
    'xdata',[],'ydata',[]);


function [flag,x1,x2] = scgedit(name,xlab1,x1,xlab2,x2)
%   Pops up a figure window to edit some data.

nvar = 1;
if nargin==5
  nvar = 2;
end

fig2 = figure('vis','off','NumberTitle','off','menu','none','name',name);
pos = get(fig2,'pos');
height = 150-34*(nvar==1);
set(fig2,'pos',[pos(1),pos(2),400,height])
uicontrol(fig2,'style','frame','pos',[0,0,400,height]);
for j=1:nvar
  height = 64+34*(nvar-j);
  eval(sprintf('xlab=xlab%d;x=x%d(:);',j,j))
  uicontrol(fig2,'style','text','pos',[10,height,80,24],...
      'string',xlab,'hor','r');
  xstr = ['[ ',sprintf(' %.4g%+.4gi ',[real(x)';imag(x)']),']'];
  ui(j) = uicontrol(fig2,'style','edit','pos',[100,height,290,24],...
      'string',xstr,'tag','scgeditbox');
end
uicontrol(fig2,'style','push','pos',[40,20,60,24],...
    'string','Clear','call',...
    'set(findobj(gcf,''tag'',''scgeditbox''),''string'',''[ ]'')')
uicontrol(fig2,'style','push','pos',[170,20,60,24],...
    'string','Done','call','global SC_FLAG,SC_FLAG=1;')
uicontrol(fig2,'style','push','pos',[300,20,60,24],...
    'string','Cancel','call','global SC_FLAG,SC_FLAG=-1;')
set(fig2,'vis','on')
global SC_FLAG
SC_FLAG = 0;
while ~SC_FLAG
  drawnow
end
flag = SC_FLAG;
clear global SC_FLAG

if flag > 0
  for j=1:nvar
    x = eval(get(ui(j),'string'),'[]');
    eval(sprintf('x%d=x;',j))
  end
end

delete(fig2)

%=======================================================
function lapsolver(varargin)

data = guidata(gcbf);
if ~isa(data.polygon,'polygon')
  errordlg('You must first define a polygon.','No polygon');
else
  p = data.polygon;
  ls = openfig('lsprogress');
  h = guihandles(ls);
  set(h.BC,'background',[1 1 0])
  bdata = lapsolvegui(p);
  set(h.BC,'background',[0 1 0])
  set(h.SC,'background',[1 1 0])
  drawnow
  phi = lapsolve(p,bdata);
  set(h.SC,'background',[0 1 0])
  set(h.Eval,'background',[1 1 0])
  drawnow
  [tri,x,y] = triangulate(p);
  u = phi(x+i*y);
  set(h.Eval,'background',[0 1 0])
  drawnow
  
  figure;
  trisurf(tri,x,y,u);
  title('Solution of Laplace''s equation')
  close(ls)
end


%=======================================================
function make_widgets(fig)

% Create uicontrol frames
set(fig,'defaultuicontrolunits','char')
% We aren't going to set positions, because that is done in resize_widgets
uicontrol('style','frame','tag','SettingsFrame');
uicontrol('style','frame','tag','ActionsFrame');


% SETTINGS
% This checkbox really seems unnecessary now.
%%setting(1) = uicontrol('Parent',fig, ...
%%    'String','Monitor progress', ...
%%    'Style','checkbox', ...
%%    'Tag','TraceBox', ...
%%    'Value',1);
setting(1) = uicontrol('Parent',fig, ...
    'String','Tolerance desired', ...
    'FontWeight','bold', ...
    'Horizontal','center',...
    'tag','TolText',...
    'Style','text');
setting(2) = uicontrol('Parent',fig, ...
    'String','1e-5', ...
    'Background',.8*[1 1 1],...
    'Style','edit', ...
    'Tag','TolEdit');
setting(3) = uicontrol('Parent',fig, ...
    'String','Canonical domain', ...
    'HorizontalAlignment','left', ...
    'Min',1, ...
    'FontWeight','bold', ...
    'horiz','cen',...
    'Style','text', ...
    'tag','DomainText',...
    'Value',1);
mapname = {'disk','half-plane','strip','rectangle','disk exterior', ...
      'disk (CR)','rectified (CR)'};
setting(4) = uicontrol('Parent',fig, ...
    'Callback',@setdomain, ...
    'Max',7, ...
    'Min',1, ...
    'String',mapname, ...
    'Style','popupmenu', ...
    'Tag','DomainPopup', ...
    'Value',1);
setting(5) = uicontrol('Parent',fig, ...
    'String','View', ...
    'HorizontalAlignment','left', ...
    'Min',1, ...
    'FontWeight','bold', ...
    'Style','text', ...
    'horiz','cen',...
    'tag','ViewText',...
    'Value',1);
setting(6) = uicontrol('Parent',fig, ...
    'Callback',@setview, ...
    'Max',3, ...
    'Min',1, ...
    'String',{'Physical domain' 'Canonical domain' 'Both domains'}, ...
    'Style','popupmenu', ...
    'Tag','ViewPopup', ...
    'Value',1);
setting(7) = uicontrol('Parent',fig, ...
    'HorizontalAlignment','right', ...
    'String','theta =', ...
    'Style','text', ...
    'horiz','left',...
    'Tag','YTText');
setting(8) = uicontrol('Parent',fig, ...
    'String','defaults', ...
    'Background',.8*[1 1 1],...
    'Style','edit', ...
    'horiz','left',...
    'Tag','YTEdit');
setting(9) = uicontrol('Parent',fig, ...
    'HorizontalAlignment','right', ...
    'String','r =', ...
    'Style','text', ...
    'horiz','left',...
    'Tag','XRText');
setting(10) = uicontrol('Parent',fig, ...
    'String','defaults', ...
    'Style','edit', ...
    'horiz','left',...
    'Background',.8*[1 1 1],...
    'Tag','XREdit');
setting(11) = uicontrol('Parent',fig, ...
    'HorizontalAlignment','left', ...
    'FontWeight','bold', ...
    'String','Mesh lines', ...
    'Style','text',...
    'horiz','cen',...
    'tag','MeshText');
 
% ACTIONS
% When loading GUI icons, we will use the system's default BG color.
bg = get(0,'defaultuicontrolbackground');
names = {'pencil','pencileraser','solve','bullseye','plot','magglass',...
      'prompt','quit'};
for j=1:length(names)
  icons{j} = imread(['private/' names{j} '.png'],'back',bg);
end
action(1) = uicontrol('parent',fig,'style','push',...
    'cdata',icons{1},...
    'tooltip','Draw new polygon',...
    'unit','pix','pos',[0 0 42 42],...
    'Callback',@draw,...
    'tag','DrawButton');
% Determine how large the icons are in character units on this computer. 
set(action(1),'unit','char')
icon_size = get(action(1),'pos');
icon_size = icon_size(3:4);
action(2) = uicontrol('parent',fig,'style','push',...
    'cdata',icons{2},...
    'tooltip','Modify current polygon',...
    'Callback',@domodify,...
    'tag','ModifyButton');
action(3) = uicontrol('parent',fig,'style','push',...
    'cdata',icons{3},...
    'tooltip','Solve for parameters',...
    'Callback',@solve,...
    'tag','SolveButton');
action(4) = uicontrol('parent',fig,'style','push',...
    'cdata',icons{4},...
    'tooltip','Set conformal center',...
    'Callback',@recenter,...
    'tag','CenterButton');
action(5) = uicontrol('parent',fig,'style','push',...
    'cdata',icons{5},...
    'tooltip','Plot a mesh',...
    'Callback',@meshplot,...
    'tag','PlotButton');
action(6) = uicontrol('parent',fig,'style','push',...
    'cdata',icons{6},...
    'tooltip','Inspect map details',...
    'Callback',@scdisplay,...
    'tag','ShowButton');
action(7) = uicontrol('parent',fig,'style','push',...
    'cdata',icons{7},...
    'tooltip','Import/export to command line',...
    'Callback',@importexport,...
    'tag','ImportExportButton');
action(8) = uicontrol('parent',fig,'style','push',...
    'cdata',icons{8},...
    'tooltip','Quit SC GUI',...
    'Callback',@scquit,...
    'tag','QuitButton');

% CONTEXT MENUS
cm = uicontextmenu('tag','PolygonContextMenu');
uimenu(cm,'label','Modify...','call',@domodify)
uimenu(cm,'label','Inspect vertices...','call',@edit)
uimenu(cm,'label','Delete...','call',@menudeletepoly)
uimenu(cm,'label','Cancel action','separator','on')

cm = uicontextmenu('tag','PointContextMenu');
uimenu(cm,'label','Inspect...','call',@editpts)
uimenu(cm,'label','Delete...','call',@deletepoints)
uimenu(cm,'label','Cancel action','separator','on')

% FIGURE MENUS
filemenu = uimenu(fig,'label','File');
uimenu(filemenu,'label','Save map data...','call',@savedata,...
       'tag','FileSaveMenu');
uimenu(filemenu,'label','Load map data...','call',@loaddata,...
       'tag','FileLoadMenu');
uimenu(filemenu,'label','Import/export to workspace...',...
       'call',@importexport,'tag','FileImportExportMenu');
uimenu(filemenu,'label','Make a printable copy','call',@printablecopy);
uimenu(filemenu,'label','Quit SC mapping...','call',@scquit);
polymenu = uimenu(fig,'label','Polygon');
uimenu(polymenu,'label','New...','call',@draw,'tag','PolyNewMenu');
uimenu(polymenu,'label','Modify...','call',@domodify,'tag','PolyModifyMenu');
uimenu(polymenu,'label','Edit numerically...','call',@edit,...
       'tag','PolyEditMenu');
mapmenu = uimenu(fig,'label','Map');
uimenu(mapmenu,'label','Solve for parameters','call',@solve,...
       'tag','MapSolveMenu');
uimenu(mapmenu,'label','Reset conformal center...','call',@recenter,...
       'tag','MapCenterMenu');
uimenu(mapmenu,'label','Plot a mesh','call',@meshplot,...
       'tag','MapPlotMenu');
uimenu(mapmenu,'label','Inspect details','call',@scdisplay,...
       'tag','MapShowMenu');
toolsmenu = uimenu(fig,'label','Tools');
tutmenu = uimenu(toolsmenu,'label','Tutorials','separator','on');
uimenu(tutmenu,'label','Basic','call','scdtutor');
uimenu(tutmenu,'label','Infinite vertices','call','scdinf');
uimenu(tutmenu,'label','Elongated polygons','call','scdlong');
uimenu(tutmenu,'label','Faber polynomials','call','scdfaber');
uimenu(toolsmenu,'label','Laplace solver...','separator','on',...
       'call',@lapsolver,'tag','ToolsLaplaceMenu');

data = guihandles(fig);

% Store all needed data
data.polygon = [];
data.map = [];
data.iscurrent = 0;
data.mapclass = { 'diskmap','hplmap','stripmap','rectmap','extermap', ...
    'crdiskmap','crrectmap' };
data.phypoints = [];
data.canpoints = [];
data.PhysicalPoints = [];
data.CanonicalPoints = [];
data.iconsize = icon_size;
data.settings = setting;
data.actions = action;
data.SCfig = fig;

% Find objects whose validity depends on polygon or map existence.
polydep = {'ModifyButton','SolveButton','FileSaveMenu','PolyModifyMenu',...
           'PolyEditMenu','MapSolveMenu','ToolsLaplaceMenu'};
mapdep = {'CenterButton','PlotButton','ShowButton','MapCenterMenu',...
          'MapPlotMenu','MapShowMenu'};
for k = 1:length(polydep)
  PD(k) = findobj(fig,'tag',polydep{k});
end
for k = 1:length(mapdep)
  MD(k) = findobj(fig,'tag',mapdep{k});
end
data.PolygonEnabled = PD(:);
data.MapEnabled = MD(:);
set(PD,'enable','off'), set(MD,'enable','off')

guidata(fig,data)

%=======================================================
function resize_widgets(fig)
% This function is called whenever the figure is resized.

figpos = get(fig,'pos');
data = guidata(fig);

% Minimum sizes
figpos(3) = max(figpos(3),95);
figpos(4) = max(figpos(4),36);

set(fig,'pos',figpos)
%movegui(fig)
h = figpos(4) - data.iconsize(2) - 0.5;
set(data.SettingsFrame,'pos',[figpos(3)-23 0 23 h])
set(data.ActionsFrame,'pos',[0 h figpos(3) data.iconsize(2)+0.5]);

% Leave a cushion around the axes
h = 10 + data.iconsize(2);
data.axeswindow = [9 5 figpos(3)-41 figpos(4)-h];

% SETTINGS
h = figpos(4) - 28 - data.iconsize(2);
w = figpos(3) - 22;
set(data.TolText, 'pos',[w+1 h+23.5 20 1.25])
set(data.TolEdit, 'pos',[w+6.5 h+21.5 9 1.5])
%set(data.TraceBox, 'pos',[w+1 h+18.75 20 1.75])
set(data.DomainText, 'pos',[w+1 h+15.5 20 1.25])
set(data.DomainPopup, 'pos',[w+1 h+13.25 20 1.75])
set(data.ViewText, 'pos',[w+1 h+10 20 1.25])
set(data.ViewPopup, 'pos',[w+1 h+7.75 20 1.75])
set(data.MeshText, 'pos',[w+1 h+4.5 20 1.25])
set(data.XRText, 'pos',[w+1 h+2.75 7 1.25])
set(data.XREdit, 'pos',[w+1 h+1 20 1.5])
set(data.YTText, 'pos',[w+1 h-1.25 7 1.25])
set(data.YTEdit, 'pos',[w+1 h-3 20 1.5])

% ACTIONS
h = figpos(4) - data.iconsize(2) - 0.25;
dx = data.iconsize(1);
dy = data.iconsize(2);
set(data.DrawButton, 'pos',[1 h dx dy])
set(data.ModifyButton, 'pos',[1+dx h dx dy])
set(data.SolveButton, 'pos',[4+2*dx h dx dy])
set(data.CenterButton, 'pos',[4+3*dx h dx dy])
set(data.PlotButton, 'pos',[7+4*dx h dx dy])
set(data.ShowButton, 'pos',[7+5*dx h dx dy])
set(data.ImportExportButton, 'pos',[10+6*dx h dx dy])
set(data.QuitButton, 'pos',[13+7*dx h dx dy])

guidata(fig,data)
% This makes an extra "flash" at startup. But it might be necessary
setview(fig)

%============================================================
function fig = make_importexportfig()

fig = figure('vis','off','NumberTitle','off','name','Import/Export',...
    'tag','sc_importexport','handlevisibility','callback',...
    'unit','char','menu','none','integerhandle','off','resize','off');

set(fig,'pos',[0,0,41,10])
movegui(fig,'center');

set(fig,'defaultuicontrolunits','char')

uicontrol(fig,'style','frame','pos',[0,0,41,10])

uicontrol('Parent',fig, ...
	'HorizontalAlignment','left', ...
	'Position',[4 6.75 18 1.5],...
	'String','Polygon variable: ', ...
	'Style','text');
uicontrol('Parent',fig, ...
	'Position',[21 6.75 10 1.5],...
	'String','', ...
	'Style','edit', ...
	'Tag','PolygonEdit');
uicontrol('Parent',fig, ...
	'HorizontalAlignment','left', ...
	'Position',[4 4.25 18 1.5],...
	'String','Map variable: ', ...
	'Style','text');
uicontrol('Parent',fig, ...
	'Position',[21 4.25 10 1.5],...
	'String','', ...
	'Style','edit', ...
	'Tag','MapEdit');

% These buttons cause the actions to be taken
%cb = ['fig=gcf; scfig = getuprop(fig,''scfig'');',...
%      'data = scgimprt(getuprop(scfig,''mapdata''));,',...
%      'setuprop(scfig,''mapdata'',data),',...
%      'scgui(scfig,''importupdate'');',...
%      'uiresume(fig)'];
uicontrol('Parent',fig, ...
	'Callback',@scimport, ...
	'Position',[3 1 9 1.75],...
	'String','Import');
%cb = ['scfig = getuprop(gcf,''scfig'');',...
%      'scgexprt(getuprop(scfig,''mapdata'')),',...
%      'uiresume(gcf)'];
uicontrol('Parent',fig, ...
	'Callback', @scexport, ...
	'Position',[16 1 9 1.75],...
	'String','Export');
uicontrol('Parent',fig, ...
	'Callback','uiresume(gcf)', ...
	'Position',[29 1 9 1.75],...
	'String','Cancel');

data = guihandles(fig);
guidata(fig,data)