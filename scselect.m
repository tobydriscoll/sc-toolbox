function K = scselect(w,beta,m,titl,msg)
%SCSELECT Select one or more vertices in a polygon.
%   K = SCSELECT(W,BETA,M) draws the polygon given by W and BETA into
%   the current figure window and then allows the user to select M
%   vertices using the mouse.  If M is not given, it defaults to 1.  On
%   exit K is a vector of indices into W.
%   
%   SCSELECT(W,BETA,M,TITLE,MESSAGE) pops up a new figure window for the
%   selection, with title TITLE and instructional message MESSAGE.
%
%   See also DRAWPOLY, PLOTPOLY, MODPOLY.

%   Copyright 1998--2001 by Toby Driscoll.
%   $Id: scselect.m 298 2009-09-15 14:36:37Z driscoll $

n = length(w);
if any(isinf(w) & isinf(w([2:n,1])))
  error('Infinite vertices must not be adjacent')
end

if nargin > 3
  fig = figure('name',titl,'numbertitle','off','integerhandle','off',...
      'menubar','none','unit','char');
  movegui(fig,'center')
  figpos = get(fig,'pos');
  pos = [0 figpos(4)-3 figpos(3) 2.5];
  t = uicontrol('style','text','unit','char','pos',pos);
  if ~iscell(msg), msg = {msg}; end
  [str,pos2] = textwrap(t,msg);
  pos(3) = pos2(3);
  pos(1) = max(0, (figpos(3)-pos2(3))/2 );
  set(t,'pos',pos,'string',str)
  h = figpos(4)-11;
  axes('unit','char','pos',[0 4 figpos(3) h]);
else
  fig = gcf;
end

[ehan,lhan] = plotpoly(w,beta,1);
turn_off_hold = ~ishold;
hold on

h = lhan;
colors = get(gca,'colororder');
if colors(1,1) > colors(1,3)
  hilit = [0 0 1];
else
  hilit = [1 0 0];
end

oldptr = get(gcf,'pointer');
set(gcf,'pointer','circle');

if nargin < 3
  m = 1;
end

% Begin selection(s)
figure(fig)
for j = 1:m
  k = [];
  while isempty(k)
    waitforbuttonpress;
    obj = get(gcf,'currentobj');
    [k,tmp] = find(obj==h);
    if isempty(k)
      disp('Selected object not a vertex.  Try again.')
    end
  end
  set(h(k,:),'color',hilit)
  drawnow
  K(j) = k;
end
set(gcf,'pointer',oldptr)

% Clean up
delete(h)
drawnow
if turn_off_hold
  hold off
end 

if nargin > 3
  delete(fig)
end
