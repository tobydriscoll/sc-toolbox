function p = polyedit(varargin)
%POLYEDIT Polygon editor.
%   P = POLYEDIT opens a window and runs the Polygon Editor. This editor
%   allows you to draw and modify a polygon using mouse input. 
%   
%   The editor has 4 main mutually exclusive modes:
%      ADD:    Add vertices in sequence by clicking the left mouse
%              button. If you hold the button down, a dashed line will
%              preview the new edge until you release. The last vertex
%              in sequence is always colored red for reference. Infinite
%              vertices may be added by clicking twice outside the axes
%              box (to indicate exit and reentry directions).
%      MOVE:   Move a vertex by clicking and dragging to a new
%              location. Only finite vertices may be moved, and they
%              must stay within the axes box. 
%      INSERT: Insert a trivial vertex by clicking on an existing
%              side. You may insert onto any side, including an infinite
%              one. 
%      DELETE: Delete a vertex by clicking on it. You may delete only
%              finite vertices. Deleting a vertex with an infinite
%              neighbor changes the angle at infinity.
%              
%   To close a polygon, draw the final side by clicking on the first
%   vertex in ADD mode or by using the appropriate pushbutton. Once a
%   polygon is closed, it remains so; adding in sequence is no longer
%   possible. You may continue to move, insert, and delete
%   vertices. When you are ready to exit, click on "OK" and the final
%   edition of your polygon will be returned in the variable P. 
%   
%   The controls on the right panel of the window allow you to constrain
%   the vertices to a grid or to discretize the allowable angles and/or
%   lengths. 
%   
%   Q = POLYEDIT(P) starts the editor with the given polygon P loaded.
%   
%   POLYEDIT(H), where H is the handle of a figure, runs the editor in
%   the given figure window. This allows you, for example, to "trace"
%   graphics you already have displayed. To do this, make sure you enter
%   HOLD ON before running POLYEDIT.
%   
%   See also the POLYGON class.

%   Copyright 2001 by Toby Driscoll.
%   $Id: polyedit.m 298 2009-09-15 14:36:37Z driscoll $

fig = [];
poly = [];

for j=1:nargin
  arg = varargin{j};
  if isgraphics(arg,'figure')  
    fig = arg;
  elseif isa(arg,'polygon')
    poly = arg;
  end
end

if isempty(fig) 
  % Open a new figure.
  fig = figure('vis','off');
  figstate = [];
  axlim = [-4 4 -4 4];
  % Make a reasonable initial size.
  un = get(0,'unit');
  set(0,'unit','char')
  ss = get(0,'screensize');
  set(0,'unit',un)
  siz = max( 0.65*ss(3:4), [95 36] );
  set(fig,'unit','char','pos',[ss(3:4)-siz-[8 4] siz])
  set(fig,'closereq',@PEQuit)
else
  % Save current state.
  turn_off_hold = ~ishold;
  figprops = {'unit';'pos';'resizefcn';'pointer';'windowbuttondown';...
	'windowbuttonup';'pointer';'menubar';'name';'numbertitle'};
  figstate = get(fig,figprops);
  axprops = {'unit';'pos';'xgrid';'ygrid'};
  axstate = get(gca,axprops);

  % Set up.
  if ~ishold
    cla
    axis auto
  end
  % Reset axes limits?
  if strcmp(get(gca,'xlimmode'),'auto') && strcmp(get(gca,'ylimmode'),'auto')
    axlim = [-4 4 -4 4];
  else
    axlim = axis;
  end
end

% Set up controls
PE_make_widgets(fig,axlim);
PE_resize_widgets(fig);
data = guidata(fig);

% Initialize with the given polygon, if any.
if ~isempty(poly)
  w = vertex(poly);
  data.polyvertex = w;
  data.polybeta = angle(poly) - 1;
  data.isclosed = 1;
  data.atinf = zeros(length(w),2);
  data.atinf(isinf(w),:) = 1;
  data.atinf = data.atinf(:);
  data.Edges = plot(poly);
  hold on
  x = get(data.Edges,'xdata');
  y = get(data.Edges,'ydata');
  for j = 1:length(w)
    data.Vertices(j) = plot(real(w(j)),imag(w(j)),'bo',...
	'markerfacecolor','b');
    if isinf(w(j))
      set(data.Vertices,'user',[ x{j} y{j} ])
    else
      set(data.Vertices,'user',[ x{j}(1) y{j}(1) ] )
    end
  end
  data.addmode = 'normal';
  guidata(fig,data)
  set(data.AddButton,'enable','off')
  set(data.MoveButton,'value',1)
  set(data.FinishButton,'enable','off')
  set(data.AcceptButton,'enable','on')
  PEMoveMode(fig)
else
  set(data.AddButton,'value',1)
  PEAddMode(fig)
end

% Let it run (modally, since it has an output value).
set(fig, 'pointer','crosshair','menubar','none',...
    'name','Polygon Editor','numbertitle','off','resizefcn',@PE_resize_widgets)
hold on
set(fig,'vis','on')
try
  uiwait(fig)
  % Retrieve the polygon and assign to output if nonempty.
  data = guidata(fig);
  P = data.polygon;
  if ~isempty(P) || nargout > 0
    p = P;
  end
catch
  errordlg({'Unexpected termination. Last message:',lasterr},...
      'PolyEdit error')
  P = [];
end
  
% Clean up the mess.
if isempty(figstate)
  % This was a new figure.
  delete(fig)
else
  set(fig,figprops,figstate)
  set(gca,axprops,axstate)
  set(data.Vertices)
  set(data.Edges)
  PEClear(fig)
  % This could be destructive in some cases.
  delete(findobj(fig,'type','uicontrol'))
  delete(data.PreviewLine)
  if turn_off_hold
    hold off
  end
  plot(p)
  set(gca,'xtickmode','auto','ytickmode','auto')
  set(gca,'xticklabelmode','auto','yticklabelmode','auto')
end


% Callback functions.

function PEAddMode(obj,varargin)
data = guidata(obj);
% Can't turn this button off by pushing on it.
if get(data.AddButton,'value')==0
  set(data.AddButton,'value',1)
  return
end
data.mode = 'Add';
guidata(obj,data)

set(data.InsertButton,'value',0)
set(data.MoveButton,'value',0)
set(data.DeleteButton,'value',0)

set(data.SnapToggle,'enable','on')
set(data.DiscAngToggle,'enable','on')
set(data.DiscLenToggle,'enable','on')
PESnapMode(obj)
PEDiscAngMode(obj)
PEDiscLenMode(obj)

set(data.Figure,'windowbuttondown',@PEAddStart)
set(data.Vertices,'button','','hittest','off')
set(data.Edges,'button','','hittest','off')


function PEInsertMode(obj,varargin)
data = guidata(obj);
% Can't turn this button off by pushing on it.
if get(data.InsertButton,'value')==0
  set(data.InsertButton,'value',1)
  return
end
data.mode = 'Insert';
guidata(obj,data)

set(data.AddButton,'value',0)
set(data.MoveButton,'value',0)
set(data.DeleteButton,'value',0)

set(data.SnapToggle,'enable','off')
set(data.DiscAngToggle,'enable','off')
set(data.DiscLenToggle,'enable','off')
PESnapMode(obj)
PEDiscAngMode(obj)
PEDiscLenMode(obj)

set(data.Figure,'windowbuttondown','')
set(data.Vertices,'button','','hittest','off')
set(data.Edges,'button',@PEInsert,'hittest','on')


function PEMoveMode(obj,varargin)
data = guidata(obj);
% Can't turn this button off by pushing on it.
if get(data.MoveButton,'value')==0
  set(data.MoveButton,'value',1)
  return
end
data.mode = 'Move';
guidata(obj,data)

set(data.AddButton,'value',0)
set(data.InsertButton,'value',0)
set(data.DeleteButton,'value',0)

set(data.SnapToggle,'enable','on')
set(data.DiscAngToggle,'enable','off')
set(data.DiscLenToggle,'enable','off')
PESnapMode(obj)
PEDiscAngMode(obj)
PEDiscLenMode(obj)

set(data.Figure,'windowbuttondown','')
set(data.Vertices,'button',@PEMoveStart,'hittest','on')
set(data.Edges,'button','','hittest','off')


function PEDeleteMode(obj,varargin)
data = guidata(obj);
% Can't turn this button off by pushing on it.
if get(data.DeleteButton,'value')==0
  set(data.DeleteButton,'value',1)
  return
end
data.mode = 'Delete';
guidata(obj,data)

set(data.AddButton,'value',0)
set(data.InsertButton,'value',0)
set(data.MoveButton,'value',0)

set(data.SnapToggle,'enable','off')
set(data.DiscAngToggle,'enable','off')
set(data.DiscLenToggle,'enable','off')
PESnapMode(obj)
PEDiscAngMode(obj)
PEDiscLenMode(obj)

set(data.Figure,'windowbuttondown','')
set(data.Vertices,'button',@PEDelete,'hittest','on')
set(data.Edges,'button','','hittest','off')


function PESnapMode(obj,varargin)
% What to do if snap mode is enabled/disabled.
data = guidata(obj);
if get(data.SnapToggle,'value') & strcmp(get(data.SnapToggle,'enable'),'on')
  % Turns off discretization modes.
  set(data.DiscAngToggle,'value',0)
  set(data.DiscLenToggle,'value',0)
  space = str2num(get(data.GridSpaceEdit,'string'));
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

function PEDiscAngMode(obj,varargin)
data = guidata(obj);
if get(data.DiscAngToggle,'value') & ...
      strcmp(get(data.DiscAngToggle,'enable'),'on')
  set(data.SnapToggle,'value',0)
  PESnapMode(obj)
end

function PEDiscLenMode(obj,varargin)
data = guidata(obj);
if get(data.DiscLenToggle,'value') & ...
      strcmp(get(data.DiscLenToggle,'enable'),'on')
  set(data.SnapToggle,'value',0)
  PESnapMode(obj)
end

function PEAddStart(obj,varargin)
data = guidata(obj);
%%set(data.PreviewLine,'vis','on')
set(data.Figure,'windowbuttonmotion',@PEAddMove)
set(data.Figure,'windowbuttonup',@PEAddFinish)
PEAddMove(obj)

function PEAddMove(obj,varargin)
data = guidata(obj);
ptrpos = get(data.Axes,'currentpoint');
P = ptrpos(1,1:2);
axlim = axis;
pts = PE_get_clicks(data);
%%x = get(data.Vertices,'xdata');
%%y = get(data.Vertices,'ydata');
%%pts = [x(:),y(:)];
m = size(pts,1);

grid = str2num(get(data.GridSpaceEdit,'string'));
grid = grid * get(data.SnapToggle,'value');
qang = 1/str2num(get(data.AngQuantEdit,'string'));
qang = qang * get(data.DiscAngToggle,'value');
qlen = str2num(get(data.LenQuantEdit,'string'));
qlen = qlen * get(data.DiscLenToggle,'value');

switch data.addmode
  case 'normal'
    mode = 0;
  case 'first'
    % Point may not be outside axes box.
    P(1) = min(max(P(1),axlim(1)),axlim(2));
    P(2) = min(max(P(2),axlim(3)),axlim(4));
    mode = 0;
  case 'infinite'
    % Point may not be inside axes box.
    if all(P>axlim([1,3])) & all(P<axlim([2,4]))
      [junk,j] = min(abs([P(1)-axlim(1:2);P(2)-axlim(3:4)]'));
      P = axlim(j);
    end
    mode = 1;
  case 'return'
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
if m > 0 && (mode~=1)
  %clearpoints(data.Previewline)
  set(data.PreviewLine,'xdata',[pts(m,1),P(1)],'ydata',[pts(m,2),P(2)]);
end
data.currentpoint = P;
guidata(obj,data)
drawnow
  
function PEAddFinish(obj,varargin)
data = guidata(obj);
P = data.currentpoint;
if ~isnan(P(1))
  set(data.Figure,'windowbuttonmotion','')
  set(data.Figure,'windowbuttonup','')
  %clearpoints(data.Previewline)
  set(data.PreviewLine,'xdata',P(1),'ydata',P(2))
  pts = PE_get_clicks(data);
  x = pts(:,1);  y = pts(:,2);
  %%x = get(data.Vertices,'xdata');
  %%y = get(data.Vertices,'ydata');
  m = length(x);
  n = m - sum(data.atinf)/2;
  axlim = axis;
  reflen = mean([diff(axlim(1:2)),diff(axlim(3:4))])/60;
  if (m > 0) && norm([P(1)-x(1), P(2)-y(1)]) < reflen
    %%x = x(1:m);
    %%y = y(1:m);
    %%set(data.Vertices,'xdata',x,'ydata',y)
    guidata(obj,data)
    PEFinish(obj)
    return
  end
  x(m+1,1) = P(1);  
  y(m+1,1) = P(2);
  
  % Beta (angle) computation.
  if n==1
    data.orient = sign( diff( x(1:2)+i*y(1:2) ) );
  elseif n > 1 
    neworient = sign( diff( x(m:m+1)+1i*y(m:m+1) ) );
    if strcmp(data.addmode,'normal')
      % Last vertex was finite.
      data.polybeta(n) = angle(data.orient/neworient) / pi;
    elseif strcmp(data.addmode,'return')
      % Last vertex was infinite.
      b = scangle( x(m-2:m+1)+1i*y(m-2:m+1) );
      b = b(2:3);
      b(b>0) = b(b>0) - 2;
      data.polybeta(n) = sum(b);
    end
    data.orient = neworient;
  end
 
  switch data.addmode
    case {'normal','first'}
      if P(1)>=axlim(1) && P(1)<=axlim(2) && P(2)>=axlim(3) && P(2)<=axlim(4)
        % Finite vertex.
        data.atinf(m+1) = 0;
        data.addmode = 'normal';
        data.polyvertex(n+1) = P(1)+1i*P(2);
        data.polybeta(n+1) = NaN;
        data.Vertices(n+1) = plot(P(1),P(2),'o','user',[P(1) P(2)]);
      else
	    data.atinf(m+1) = 1;
        data.addmode = 'infinite';
        % Set up for re-entry point next time.
        set(data.Figure,'pointer','cross')
        h = [data.InsertButton,data.MoveButton,data.DeleteButton,...
            data.FinishButton];
        set(h,'enable','off')
        data.polyvertex(n+1) = Inf;
        data.polybeta(n+1) = NaN;
        data.Vertices(n+1) = plot(Inf,Inf,'o','user',[P(1) P(2)]);
      end
      if n > 0
        data.Edges(n) = plot(x(m:m+1),y(m:m+1),'-');
        set(data.Edges(n),'user',n)
      end
    case 'infinite'
      data.atinf(m+1) = 1;
      set(data.Figure,'pointer','crosshair')
      data.addmode = 'return';
      pt = get(data.Vertices(n+0.5),'user');
      set(data.Vertices(n+0.5),'user',[pt;[P(1) P(2)]])
    case 'return'
      data.Edges(n) = plot(x(m:m+1),y(m:m+1),'-');
      set(data.Edges(n),'user',n)
      data.atinf(m+1) = 0;
      data.polyvertex(n+1) = P(1)+1i*P(2);
      data.polybeta(n+1) = NaN;
      data.Vertices(n+1) = plot(P(1),P(2),'o','user',[P(1) P(2)]);
      % Restore mode buttons.
      h = [data.InsertButton,data.MoveButton,data.DeleteButton,...
	    data.FinishButton];
      set(h,'enable','on')
      data.addmode = 'normal';
  end     
  
  set(data.Vertices(1:end-1),'markerfacecolor','b','color','b')
  set(data.Vertices(end),'markerfacecolor','r','color','r')
  data.currentpoint = [NaN NaN];
  guidata(obj,data)
  drawnow
  
end

function PEMoveStart(obj,varargin)
data = guidata(obj);
% Figure out which vertex is selected.
%%P = get(gca,'currentpoint');
%%[m,idx] = min( abs( P(1)+i*P(2) - data.polyvertex ) );
%%data.moveselected = idx(1);
idx = find( get(gcf,'currentobj')==data.Vertices );
data.moveselected = idx;
adj = PE_adjacent_edges(data,data.moveselected);
data.moveadj = adj;
guidata(obj,data)
set(data.Edges(adj(~isnan(adj))),'linest','--')
set(data.Figure,'windowbuttonmotion',@PEMoveMove)
set(data.Figure,'windowbuttonup',@PEMoveFinish)

function PEMoveMove(obj,varargin)
data = guidata(obj);
xy = get(gca,'currentpoint');
xy = xy(1,1:2);
w = data.polyvertex;
n = length(w);
% Keep inside axes box.
axlim = axis;
xy(1) = min( max(xy(1),axlim(1)), axlim(2) );
xy(2) = min( max(xy(2),axlim(3)), axlim(4) );
% Snap to grid if requested.
if get(data.SnapToggle,'value');
  grid = str2num(get(data.GridSpaceEdit,'string'))*[1 1];
  minxy = axlim([1 3]);
  xy = minxy + grid.*(round((xy-minxy)./grid));
end
k = data.moveselected;
set(data.Vertices(k),'xdata',xy(1),'ydata',xy(2),'user',xy)

adj = data.moveadj;
edges = data.Edges;
beta = data.polybeta;

scale = max( axlim(2)-axlim(1), axlim(4)-axlim(3) );
if ~isnan(adj(1))
  % Adjust the "predecessor" side.
  j = rem(k-2+n,n)+1;
  if isinf(w(j))
    xd = get(edges(j),'xd');
    yd = get(edges(j),'yd');
    phi = atan2(-diff(yd),-diff(xd));
    % Note: draw lines to infinity long enough
    r = 2*scale;
    %r = sqrt(diff(xd)^2+diff(yd)^2);
    y = xy(1)+1i*xy(2) + [r*exp(1i*phi),0];
    % The "clickpoint" for the inf vertex has changed.
    pts = get(data.Vertices(j),'user');
    pts(2,:) = [real(y(1)) imag(y(1))];
    set(data.Vertices(j),'user',pts)
  else
    y = [w(j),xy(1)+1i*xy(2)];
    phi = angle(diff(y)/diff(w([j,k])));
    beta(k) = beta(k)+phi/pi;
    beta(j) = beta(j)-phi/pi;
  end
  %clearpoints(edges(adj(1)))
  set(edges(adj(1)),'xdata',real(y),'ydata',imag(y))
end

if ~isnan(adj(2))
  % Adjust the "successor" side.
  j = rem(k,n)+1;
  if isinf(w(j))
    [xd,yd] = getpoints(edges(k));
    phi = atan2(diff(yd),diff(xd));
    r = 2*scale;
    %%r = sqrt(diff(xd)^2+diff(yd)^2);
    y = xy(1)+1i*xy(2) + [0,r*exp(1i*phi)];
    % The "clickpoint" for the inf vertex has changed.
    pts = get(data.Vertices(j),'user');
    pts(1,:) = [real(y(2)) imag(y(2))];
    set(data.Vertices(j),'user',pts)
  else
    y = [xy(1)+1i*xy(2),w(j)];
    phi = angle(-diff(y)/diff(w([j,k])));
    beta(k) = beta(k)-phi/pi;
    beta(j) = beta(j)+phi/pi;
  end
  set(edges(adj(2)),'xdata',real(y),'ydata',imag(y))
end
% Make changes effective
drawnow
data.polyvertex(k) = xy(1)+1i*xy(2);
data.polybeta = beta;
guidata(obj,data)

function PEMoveFinish(obj,varargin)
data = guidata(obj);
set(data.Figure,'windowbuttonmotion','')
set(data.Figure,'windowbuttonup','')
set(data.Edges,'linesty','-')
% Update the orientation of the last side, if polygon is not closed.
if ~data.isclosed
  pts = PE_get_clicks(data);
  w = pts(:,1) + 1i*pts(:,2);
  data.orient = sign( diff( w(end-1:end) ) );
  guidata(obj,data)
end
    

function PEInsert(obj,varargin)
data = guidata(obj);
idx = find( get(gcf,'currentobj')==data.Edges );
[wn,bn] = scaddvtx(data.polyvertex,data.polybeta,idx,axis);
P = wn(idx+1);
data.polyvertex = wn;
data.polybeta = bn;
%%new = zeros(1,length(data.atinf)+1);
%%new([1:idx idx+2:end]) = data.atinf(:);
%%data.atinf = new;
data.atinf = [data.atinf(1:idx); 0; data.atinf(idx+1:end)];
newVertex = plot(real(P),imag(P),'bo',...
    'markerfacecolor','b','user',[real(P) imag(P)]);
data.Vertices = [data.Vertices(1:idx) newVertex data.Vertices(idx+1:end)];
pred = data.Edges(idx);
x = get(pred,'xdata');  y = get(pred,'ydata');
set(pred,'xdata',[x(1) real(P)],'ydata',[y(1) imag(P)])
succ = plot([real(P) x(2)],[imag(P) y(2)],'-');
data.Edges = [data.Edges(1:idx) succ data.Edges(idx+1:end)];
guidata(obj,data)


function PEDelete(obj,varargin)
data = guidata(obj);
k = find( get(gcf,'currentobj')==data.Vertices );
n = length(data.polyvertex);
if n==1
  PEClear(obj)
  return
end
% Special case for an unfinished polygon at an endpoint.
if ~data.isclosed && ( k==1 || k==n ) 
  if (k==1 && isinf(data.polyvertex(2))) || ...
	(k==n && isinf(data.polyvertex(n-1))) 
    errordlg('To be deleted, a vertex must have at least one finite neighbor.','PolyEdit Error')
    return
  end
  delete(data.Vertices(k))
  data.Vertices(k) = [];
  data.polyvertex(k) = [];
  data.polybeta(k) = [];
  data.atinf(k) = [];
  if k==1
    data.polybeta(1) = NaN;
    delete(data.Edges(1))
    data.Edges(1) = [];
  else
    data.orient = data.orient*exp(1i*pi*data.polybeta(end));
    data.polybeta(end) = NaN;
    delete(data.Edges(n-1))
    data.Edges(n-1) = [];
    set(data.Vertices(n-1),'color','r','markerfacecolor','r')
  end
else
  w = data.polyvertex;
  beta = data.polybeta;
  % Neighbors
  pred = rem( k-2+n, n ) + 1;
  succ = rem(k,n) + 1;
  % There must be at least one finite neighbor.
  infp = isinf(w(pred));
  infs = isinf(w(succ));
  if ~infp && ~infs
    wp = w(pred);
    ws = w(succ);
    beta(pred) = angle( exp(1i*pi*beta(pred))*(w(k)-wp)/(ws-wp) ) / pi;
    beta(succ) = angle( exp(1i*pi*beta(succ))*(wp-ws)/(w(k)-ws) ) / pi;
    side = [wp ws];
  elseif infp && ~infs
    axlim = axis;
    scale = 2*max( axlim(2)-axlim(1), axlim(4)-axlim(3) );
    beta(pred) = beta(pred) + beta(k);
    side = w(succ) + [scale*sign(w(k)-w(succ)) 0];
    pts = get(data.Vertices(pred),'user');
    pts(2,:) = [real(side(1)) imag(side(1))];
    set(data.Vertices(pred),'user',pts)
  elseif ~infp && infs
    axlim = axis;
    scale = 2*max( axlim(2)-axlim(1), axlim(4)-axlim(3) );
    beta(succ) = beta(succ) + beta(k);
    side = w(pred) + [0 scale*sign(w(k)-w(pred))];
    pts = get(data.Vertices(succ),'user');
    pts(1,:) = [real(side(2)) imag(side(2))];
    set(data.Vertices(succ),'user',pts)
  else
    errordlg('To be deleted, a vertex must have at least one finite neighbor.','PolyEdit Error')
    return
  end

  data.polyvertex(k) = [];
  data.atinf(k) = [];
  data.polybeta = beta;
  data.polybeta(k) = [];
  delete(data.Vertices(k))
  data.Vertices(k) = [];
  %clearpoints(data.Edges(pred))
  set(data.Edges(pred),'xdata',real(side),'ydata',imag(side))
  delete(data.Edges(k))
  data.Edges(k) = [];
end
guidata(data.Figure,data)


function PEClear(obj,varargin)
data = guidata(obj);
delete(data.Vertices)
delete(data.Edges)
data.Vertices = [];
data.Edges = [];
data.polyvertex = [];
data.polybeta = [];
data.atinf = [];
data.isclosed = 0;
data.addmode = 'first';
guidata(obj,data)
set(data.AddButton,'value',1)
PEAddMode(obj)


function PEFinish(obj,varargin)
data = guidata(obj);
pts = PE_get_clicks(data);
x = pts(:,1);  y = pts(:,2);
%%x = get(data.Vertices,'xdata');
%%y = get(data.Vertices,'ydata');
m = length(x);
n = length(x) - sum(data.atinf)/2;
if n > 1
  data.Edges(n) = plot(x([m,1]),y([m,1]),'-','user',n);
  if n==2
    % Special degenerate case.
    data.polybeta = [1;1];
  else
    % Find betas at last and first vertex.
    neworient = sign( diff( x([m 1])+i*y([m 1]) ) );
    if ~data.atinf(m)
      % Last vertex was finite.
      data.polybeta(n) = angle(data.orient/neworient) / pi;
    else
      % Last vertex was infinite.
      b = scangle( x([m-2:m 1])+i*y([m-2:m 1]) );
      b = b(2:3);
      b(b>0) = b(b>0) - 2;
      data.polybeta(n) = sum(b);
    end
    data.orient = neworient;
    % First vertex is always finite.
    neworient = sign( diff( x(1:2)+i*y(1:2) ) );
    data.polybeta(1) = angle(data.orient/neworient) / pi; 
  end
end
drawnow
% Let the user know if this thing is not legal.
try 
  polygon(data.polyvertex,data.polybeta+1);
catch
  warndlg('This is currently not a proper polygon','PolyEdit Warning',...
      'modal')
end
data.isclosed = 1;
set(data.Vertices,'color','b','markerfacecolor','b')
guidata(obj,data)
set(data.AddButton,'enable','off')
set(data.MoveButton,'value',1)
PEMoveMode(obj)
set(data.FinishButton,'enable','off')
set(data.AcceptButton,'enable','on')

function PEXLimSet(obj,varargin)
data = guidata(obj);
xlim = str2num(get(data.XLimEdit,'string'));
set(gca,'xlim',xlim)
PESnapMode(obj)
 
function PEYLimSet(obj,varargin)
data = guidata(obj);
ylim = str2num(get(data.YLimEdit,'string'));
set(gca,'ylim',ylim)
PESnapMode(obj)
 
function PEAccept(obj,varargin)
data = guidata(obj);
try 
  p = polygon( data.polyvertex, 1+data.polybeta);
catch
  errordlg('Not a proper polygon.','Polyedit error')
  p = polygon([]);
end
data.polygon = p;
guidata(obj,data)
uiresume

function PEQuit(obj,varargin)
data = guidata(obj);
p = polygon([]);
data.polygon = p;
guidata(obj,data)
uiresume


% Utility functions.

function pts = PE_get_clicks(data)
pts = get(data.Vertices,'user');
if isempty(pts)
  pts = zeros(0,2);
elseif size(pts,1)>1
  pts = cat(1,pts{:});
end

function adj = PE_adjacent_edges(data,idx)
n = length(data.polyvertex);
m = length(data.Edges);
if idx==1
  if ~data.isclosed
    adj = [NaN 1];
  else
    adj = [m 1];
  end
elseif idx==n & ~data.isclosed
  adj = [m NaN];
else
  adj = [idx-1 idx];
end


% GUI functions.

function PE_make_widgets(fig,axlim)

set(fig,'units','char','defaultuicontrolunits','char')
set(gca,'xlim',axlim(1:2),'ylim',axlim(3:4))

% Create uicontrol frames
uicontrol('style','frame','tag','SettingsFrame');
uicontrol('style','frame','tag','ActionsFrame');

% Load in icon data for settings
bg = get(0,'defaultuicontrolbackground');
names = {'addseq','addexist','move','eraser',...
      'trash','polyclose','OK','quit'};
for j=1:length(names)
  icons{j} = imread(['private/' names{j} '.png'],'back',bg);
end

widget(1) = uicontrol('style','toggle',...
    'cdata',icons{1},...
    'tooltip','Add vertices in sequence',...
    'unit','pix',...
    'pos',[0 0 42 42],...
    'call',@PEAddMode,...
    'tag','AddButton');
widget(end+1) = uicontrol('style','toggle',...
    'cdata',icons{2},...
    'tooltip','Insert vertices',...
    'unit','pix',...
    'pos',[0 0 42 42],...
    'call',@PEInsertMode,...
    'tag','InsertButton');
widget(end+1) = uicontrol('style','toggle',...
    'cdata',icons{3},...
    'tooltip','Move vertices',...
    'unit','pix',...
    'pos',[0 0 42 42],...
    'call',@PEMoveMode,...
    'tag','MoveButton');
widget(end+1) = uicontrol('style','toggle',...
    'cdata',icons{4},...
    'tooltip','Delete vertices',...
    'unit','pix',...
    'pos',[0 0 42 42],...
    'call',@PEDeleteMode,...
    'tag','DeleteButton');
widget(end+1) = uicontrol('style','push',...
    'cdata',icons{5},...
    'tooltip','Start over',...
    'unit','pix',...
    'pos',[0 0 42 42],...
    'call',@PEClear,...
    'tag','ClearButton');
widget(end+1) = uicontrol('style','push',...
    'cdata',icons{6},...
    'tooltip','Add final side',...
    'unit','pix',...
    'pos',[0 0 42 42],...
    'call',@PEFinish,...
    'tag','FinishButton');
widget(end+1) = uicontrol('style','push',...
    'cdata',icons{7},...
    'tooltip','Accept and exit',...
    'unit','pix',...
    'pos',[0 0 42 42],...
    'call',@PEAccept,...
    'enable','off',...
    'tag','AcceptButton');
widget(end+1) = uicontrol('style','push',...
    'cdata',icons{8},...
    'tooltip','Discard and quit',...
    'unit','pix',...
    'pos',[0 0 42 42],...
    'call',@PEQuit,...
    'tag','CancelButton');

names = {'snapxy','quantlen','quantang'};
for j=1:length(names)
  icons{j} = imread(['private/' names{j} '.png'],'back',bg);
end

widget(end+1) = uicontrol('style','text',...
    'string','x limits:',...
    'hor','r',...
    'tag','XLimText');
widget(end+1) = uicontrol('style','text',...
    'string','y limits:',...
    'hor','r',...
    'tag','YLimText');
widget(end+1) = uicontrol('style','edit',...
    'string',sprintf('[%.3g %.3g]',axlim(1:2)),...
    'tag','XLimEdit',...
    'call',@PEXLimSet);
widget(end+1) = uicontrol('style','edit',...
    'string',sprintf('[%.3g %.3g]',axlim(3:4)),...
    'tag','YLimEdit',...
    'call',@PEYLimSet );

widget(end+1) = uicontrol('style','toggle',...
    'cdata',icons{1},...
    'tooltip','Snap to grid',...
    'unit','pix',...
    'pos',[0 0 42 42],...
    'call',@PESnapMode,...
    'tag','SnapToggle');

widget(end+1) = uicontrol('style','toggle',...
    'cdata',icons{2},...
    'tooltip','Discretize length',...
    'unit','pix',...
    'pos',[0 0 42 42],...
    'call',@PEDiscLenMode,...
    'tag','DiscLenToggle');
widget(end+1) = uicontrol('style','toggle',...
    'cdata',icons{3},...
    'tooltip','Discretize angle',...
    'unit','pix',...
    'pos',[0 0 42 42],...
    'call',@PEDiscAngMode,...
    'tag','DiscAngToggle');

widget(end+1) = uicontrol('style','text',...
    'string','Spacing:',...
    'hor','r',...
    'tag','SpacingText');
widget(end+1) = uicontrol('style','text',...
    'string','Quantum:',...
    'hor','r',...
    'tag','AngQuantText');
widget(end+1) = uicontrol('style','text',...
    'string','Quantum:',...
    'hor','r',...
    'tag','LenQuantText');

gridspace = diff(axlim(1:2))/16;
widget(end+1) = uicontrol('style','edit',...
    'string',sprintf('%.3g',gridspace),...
    'call',@PESnapMode,...
    'tag','GridSpaceEdit');
widget(end+1) = uicontrol('style','text',...
    'string','p / ',...
    'fontname','symbol',...
    'fontsize',11,...
    'hor','r',...
    'tag','PiText');
widget(end+1) = uicontrol('style','edit',...
    'string','6',...
    'tag','AngQuantEdit');
widget(end+1) = uicontrol('style','edit',...
    'string',sprintf('%.3g',2*gridspace),...
    'tag','LenQuantEdit');

set(gca,'unit','char','box','on','xgrid','off','ygrid','off',...
     'plotboxaspectratio',[1,1,1],'tag','Axes')

% Preview line
animatedline('linestyle','--','color','r','clipping','off',...
    'tag','PreviewLine');

data = guihandles(fig);

% Determine how large the icons are in character units on this computer. 
set(widget,'unit','char')
icon_size = get(widget(1),'pos');
data.iconsize = icon_size(3:4);

data.Figure = fig;
data.Edges = [];
data.Vertices = [];
data.currentpoint = [NaN NaN];
data.addmode = 'first';
data.atinf = [];
data.isclosed = 0;

guidata(fig,data)


function PE_resize_widgets(fig,varargin)

figpos = get(fig,'pos');
data = guidata(fig);

% Minimum sizes
figpos(3) = max(figpos(3),95);
figpos(4) = max(figpos(4),36);

set(fig,'pos',figpos)

h = figpos(4) - data.iconsize(2) - 0.5;
set(data.SettingsFrame,'pos',[figpos(3)-23 0 23 h])
set(data.ActionsFrame,'pos',[0 h figpos(3) data.iconsize(2)+0.5]);

% Leave a cushion around the axes
h = 10 + data.iconsize(2);
set(data.Axes,'pos',[ 8 5 figpos(3)-39 figpos(4)-h]);

% Settings
is = data.iconsize;
w = figpos(3) - 22;
wb = (21-is(1))/2;
h = figpos(4) - 2 - 2*is(2);
set(data.SnapToggle, 'pos',[w+wb h is(1) is(2)])
h = h - 1.75;
set(data.SpacingText, 'pos',[w+1 h 10 1.25])
set(data.GridSpaceEdit, 'pos',[w+11.5 h 9.5 1.5])
h = h - 2.5 - is(2);
set(data.DiscLenToggle, 'pos',[w+wb h is] )
h = h - 1.75;
set(data.LenQuantText, 'pos',[w+1 h 11 1.25])
set(data.LenQuantEdit, 'pos',[w+13 h 6 1.5])
h = h - 2.5 - is(2);
set(data.DiscAngToggle, 'pos',[w+wb h is] )
h = h - 1.75;
set(data.AngQuantText, 'pos',[w+1 h 10 1.25])
set(data.PiText, 'pos',[w+11.5 h 4 1.25])
set(data.AngQuantEdit, 'pos',[w+16 h 5 1.5])
h = h - 4;
set(data.XLimText, 'pos',[w+1 h 9 1.25])
set(data.XLimEdit, 'pos',[w+11 h 10 1.5])
h = h - 2.5;
set(data.YLimText, 'pos',[w+1 h 9 1.25])
set(data.YLimEdit, 'pos',[w+11 h 10 1.5])

% Actions
h = figpos(4) - data.iconsize(2) - 0.25;
dx = data.iconsize(1);
dy = data.iconsize(2);
set(data.AddButton, 'pos',[1 h dx dy])
set(data.MoveButton, 'pos',[1+dx h dx dy])
set(data.InsertButton, 'pos',[1+2*dx h dx dy])
set(data.DeleteButton, 'pos',[1+3*dx h dx dy])
set(data.ClearButton, 'pos',[4+4*dx h dx dy])
set(data.FinishButton, 'pos',[4+5*dx h dx dy])
set(data.AcceptButton, 'pos',[7+6*dx h dx dy])
set(data.CancelButton, 'pos',[7+7*dx h dx dy])

drawnow
