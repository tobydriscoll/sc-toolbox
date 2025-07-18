function bdata = lapsolvegui(varargin)
% LAPSOLVEGUI GUI implemtentation for lapsolvegui.fig.

%    FIG = LAPSOLVEGUI launch lapsolvegui GUI.
%    LAPSOLVEGUI('callback_name', ...) invoke the named callback.

%    Copyright 2002 by Toby Driscoll.
%    $Id: lapsolvegui.m 298 2009-09-15 14:36:37Z driscoll $

% Thanks to Paul Syred for suggestions on this GUI concept.

if isa(varargin{1},'polygon')

  fig = openfig(mfilename,'reuse');
  
  % Use system color scheme for figure:
  bg = get(0,'defaultUicontrolBackgroundColor');
  set(fig,'Color',bg);
  
  % Generate a structure of handles to pass to callbacks, and store it. 
  handles = guihandles(fig);
  
  % Set up initial data.
  p = varargin{1};  n = length(p);
  set(fig,'handlevis','on')
  axes(handles.Axes);  hold on, box on
  side = plot(p);
  plot(vertex(p),'b.','markersize',6)
  set(fig,'handlevis','callback')
  set(side,{'tag'},cellstr(num2str((1:n)')));
  set(side,'buttondown',...
           'lapsolvegui(''PolygonSide_Callback'',gcbo,[],guidata(gcbo))');
  set(side(1),'color','r')  
      
  set(handles.ListSides,'string',cellstr(num2str(zeros(n,1))),'value',1);
  set(handles.RadioDirichlet,'value',1)
  set(handles.EditValue,'string','0')
  set(handles.RadioNeumann,'value',0)
  
  % The radio buttons wouldn't work right for me.
  set(handles.RadioDirichlet,'style','check')
  set(handles.RadioNeumann,'style','check')
  
  set(handles.ButtonOK,'background',bg.^0.5)
  
  handles.PolygonSides = side;
  guidata(fig, handles);

  % Wait for callbacks to run and window to be dismissed:
  uiwait(fig);

  % Read off the BCs and exit.
  BC = get(handles.ListSides,'string');
  bdata = NaN(n,1);
  for k = 1:n
    if ~strcmp( BC{k}, 'N' )
      bdata(k) = str2num( BC{k} );
    end
  end
  
  close(fig)
  
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
  
  try
    if (nargout)
      [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
    else
      feval(varargin{:}); % FEVAL switchyard
		end
  catch
    disp(lasterr);
  end
  
end


% --------------------------------------------------------------------
function varargout = PolygonSide_Callback(h, eventdata, handles, ...
                                          varargin)
% Check for extension of selection.
if strcmp(get(gcbf,'selectiontype'),'normal')
  target = str2num( get(h,'tag') );
else
  target = get(handles.ListSides,'value');
  new = str2num( get(h,'tag') );
  target(end+1) = new;
end

selectsides(handles,target);



% --------------------------------------------------------------------
function varargout = ListSides_Callback(h, eventdata, handles, varargin)

% Deduce the side number from the selection in the listbox.
target = get(h,'value');
selectsides(handles,target);


% --------------------------------------------------------------------
function varargout = selectsides(handles, target)

% Make sure the list box is up to date.
set(handles.ListSides,'value',target);

% Re-color every polygon side (we don't know who was selected before).
BC = get(handles.ListSides,'string');
side = handles.PolygonSides;
for k = 1:length(side)
  if strcmp(BC{k},'N')
    set(side(k),'color','k')
  else
    set(side(k),'color','b')
  end
end

% Now color selected sides.
set(side(target),'color','r')
drawnow

% Update the radio buttons and value box, unless this was a selection to
% multiple sides.
if isscalar(target)
  if strcmp(BC{target},'N')
    set(handles.RadioDirichlet,'value',0)
    set(handles.EditValue,'enable','off')
    set(handles.RadioNeumann,'value',1)
  else
    set(handles.RadioDirichlet,'value',1)
    set(handles.EditValue,'enable','on','string',BC{target})
    set(handles.RadioNeumann,'value',0)
  end
end

% --------------------------------------------------------------------
function varargout = setvalue(handles, val)

target = get(handles.ListSides,'value');
BC = get(handles.ListSides,'string');
for k = 1:length(target)
  BC{target(k)} = val;
end
set(handles.ListSides,'string',BC)

 

% --------------------------------------------------------------------
function varargout = RadioDirichlet_Callback(h, eventdata, handles, varargin)

set(handles.EditValue,'enable','on')
set(handles.RadioDirichlet,'value',1)
set(handles.RadioNeumann,'value',0)
setvalue(handles,get(handles.EditValue,'string'));


% --------------------------------------------------------------------
function varargout = RadioNeumann_Callback(h, eventdata, handles, varargin)

set(handles.EditValue,'enable','off')
set(handles.RadioDirichlet,'value',0)
set(handles.RadioNeumann,'value',1)
setvalue(handles,'N');


% --------------------------------------------------------------------
function varargout = EditValue_Callback(h, eventdata, handles, varargin)

val = get(h,'string');
if isempty(str2num(val))
  % Invalid entry; reset.
  target = get(handles.ListSides,'value');
  BC = get(handles.ListSides,'string');
  set(h,'string',BC{target(1)});
else
  % Accept entry.
  val = num2str(str2num(val));  % clean it up
  set(h,'string',val)
  setvalue(handles,val)
end


% --------------------------------------------------------------------
function varargout = ButtonOK_Callback(h, eventdata, handles, varargin)

uiresume(gcbf);