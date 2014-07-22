function [wr,betar,affr] = crrect(w,beta,cr,aff,Q,betar)
%CRRECT Graphically create a rectified map.
%   A rectified map is one between a generic polygon and a polygon
%   having all angles as multiples of pi/2. Once the CR parameter
%   problem is solved, you can define a rectified polygon to map to by
%   specifying the rectified angles. (Side lengths are not free
%   variables; they are a consequence of the original polygon and the
%   rectified angles.) There is no unique correct choice of rectified
%   angles. No automatic selection of the rectified angles is currently
%   available, so CRRECT allows you to do it graphically.
%       
%   [WR,BETAR,AFFR] = CRRECT(W,BETA,CR,AFF,Q) first selects four
%   vertices of the polygon W to serve as corners of a rectangle, which
%   is the initial rectified domain. Then a side-by-side view of the
%   rectified and original polygons is created. The vertices in each
%   view are color-coded according to angle. There is a pop-up selection
%   list that allows you to pick one rectified angle/color. Any vertex
%   you click on (in either view) will be assigned that color and angle.
%       
%   A counter shows the current sum of angles. When this is -2, you have
%   defined a valid rectified map. If you click the button marked
%   Recompute, the rectified polygon will be calculated and
%   displayed. You may make further changes. If you click "Done", CRRECT
%   returns the rectified vertices and angles, and the affine constants
%   needed to compute the rectified map at specific points. You may also
%   click "Cancel" to return nothing.
%       
%   IMPORTANT NOTE: The rectified polygon may not be embeddable in the
%   plane. That is, it may overlap itself. This may be difficult to see
%   graphically, but it does not affect CRRECT in any way.
%       
%   An additional input parameter will be used as BETAR, and no
%   graphical procedure will take place.
%       
%   See also CRPARAM, CRRMAP, CRRPLOT.
       
%   Copyright 1998 by Toby Driscoll.
%   $Id: crrect.m 7 1998-05-10 04:37:19Z tad $
n = length(w);
if nargin < 6
  % betar not given: use graphical procedure
  
  % Start with a rectangle
  % Pick the 4 vertices closest to right angles
  [tmp,k] = sort(abs(beta+.5));
  k = k(1:4);
  betar = 0*beta;
  betar(k) = -.5*ones(4,1);
  wr = NaN*w;
  wr([k(1) rem(k(1),n)+1]) = [0;1];
  [affr,wr] = crdiskmap.craffine(wr,betar,cr,Q);
  %aff = craffine(w,beta,cr,Q);
  
  % Move on to graphical refinement
  fig = figure;
  
  % Dummy objects that store the data
  uicontrol('style','frame','vis','off','user',wr,...
      'tag','wr');
  uicontrol('style','frame','vis','off','user',betar,...
      'tag','betar');
  uicontrol('style','frame','vis','off','user',cr,...
      'tag','cr');
  uicontrol('style','frame','vis','off','user',Q,...
      'tag','Q');
  uicontrol('style','frame','vis','off','user',affr,...
      'tag','affr');
  
  % Set up figure 
  pos = get(fig,'pos');
  set(fig,'pos',[pos(1:2) max(pos(3),700) pos(4)],'user',0)
  
  % Rectified plane
  ax(1) = axes('pos',[.56 .29 .4 .68],'tag','crrect_rectified','box','on');
  hold on
  [hr,h] = plotpoly(wr,betar,1);
  uicontrol('style','frame','vis','off','tag','rect_sides','user',hr(:));
  set(h(:,1),'markersize',20)
  set(h(k,:),'color','m')		% the 4 selected
  set(ax(1),'user',h)
  
  % Original (physical) plane
  ax(2) = axes('pos',[.06 .29 .4 .68],'tag','crrect_original','box','on');
  hold on
  [ho,h] = plotpoly(w,beta,1);
  uicontrol('style','frame','vis','off','tag','orig_sides','user',ho(:));
  axlim = axis;
  set(h(:,1),'markersize',20)
  set(h(k,:),'color','m')		% the 4 selected
  set(ax(2),'user',h)
  p = max(pos(3),700);
  
  set(ax,'dataaspectratio',[1 1 1],'plotboxaspectratio',[1 1 1])
  
  % GUI controls
  uicontrol('style','frame','unit','pix','pos',[0 0 p 100])

  uicontrol('style','text','unit','pix','pos',[10 48 100 16],...
      'string','Set rectified','hor','cen')
  uicontrol('style','text','unit','pix','pos',[10 33 100 16],...
      'string','angle to:','hor','cen')
  
  rh(1)=uicontrol('style','radio','unit','pix','pos',[110 72 80 22],...
      'foreground',[1 0 1],'string','pi / 2','tag','angradio');
  rh(2)=uicontrol('style','radio','unit','pix','pos',[110 49 80 22],...
      'foreground',[0 0 1],'string','pi','tag','angradio');
  rh(3)=uicontrol('style','radio','unit','pix','pos',[110 26 80 22],...
      'foreground',[0 0.5 0.5],'string','3 pi / 2','tag','angradio');
  rh(4)=uicontrol('style','radio','unit','pix','pos',[110 3 80 22],...
      'foreground',[0 0 0],'string','2 pi','tag','angradio');
  set(rh,'user',rh,'call','cback(''radioexclude'')');
  set(rh(2),'value',1)
  
  uicontrol('style','text','unit','pix','pos',[400 60 80 20],...
      'string','sum(beta) =')
  uicontrol('style','text','unit','pix','pos',[410 36 60 25],...
      'string','-2.0','tag','sumbeta')
  
  uicontrol('style','push','unit','pix','pos',[p-110 68 92 25],...
      'string','Recompute','tag','recompute','call','cback(''compute'');')
  uicontrol('style','push','unit','pix','pos',[p-110 37 92 25],...
      'string','Done','call','set(gcf,''user'',1)')
  uicontrol('style','push','unit','pix','pos',[p-110 7 92 25],...
      'string','Cancel','call','set(gcf,''user'',-1)')
  
  % Let the GUI operate until done
  %set(fig,'windowbuttondown','cback(''sideselect'');')
  set(fig,'windowbuttonup','cback(''select'');')
  while(~get(fig,'user'))
    drawnow
  end
  
  % Return values and clean up
  if get(fig,'user') > 0
    % Normal exit
    wr = get(findobj(gcf,'tag','wr'),'user');
    betar = get(findobj(gcf,'tag','betar'),'user');
    affr = get(findobj(gcf,'tag','affr'),'user');
    delete(fig)
  else
    % Cancel
    delete(fig)
    error('Procedure canceled.')
  end
  
else
  % betar was supplied: just do computation
  
  aff = crdiskmap.craffine(w,beta,cr,Q);
  k = min(find(betar<0));
  wr = NaN*w;
  wr([k rem(k,n)+1]) = [0;1];
  [affr,wr] = crdiskmap.craffine(wr,betar,cr,Q);
end

% End of user-call section

function cback(cmd)
% Called from the GUI

switch cmd
case 'select'
  % An object is selected
  % Is it a vertex? Which one?
  ax = findobj(gcf,'type','axes');
  obj = get(gcf,'currentobj');
  k = find(any((obj==get(ax(1),'user'))') | any((obj==get(ax(2),'user'))'));
  if ~isempty(k) 
    % Get current angle value and assign it
    betar = get(findobj(gcf,'tag','betar'),'user');
    %%      t = findobj(gcf,'tag','betamenu');
    %%      val = get(t,'value');
    %%      betar(k) = (val-2)/2;
    t = findobj(gcf,'tag','angradio');
    t = get(t(1),'user');
    val = get(t,'value');
    val = find( cat(1,val{:}) );
    betar(k) = (val-2) / 2;
    set(findobj(gcf,'tag','betar'),'user',betar)
    % New colors
    colortable = [1 0 1;0 0 1;0 .5 .5;0 0 0];
    h = get(ax(1),'user');
    set(h(k,:),'color',colortable(val,:))
    h = get(ax(2),'user');
    set(h(k,:),'color',colortable(val,:))
    set(findobj(gcf,'tag','sumbeta'),'string',sprintf('%.1f',sum(betar)))
    hr = findobj(gcf,'tag','recompute');
    if abs(sum(betar)+2) > 1e-10
      set(hr,'enable','off')
    else
      set(hr,'enable','on')
    end
  end
  
case 'sideselect'
  % If a side is selected, highlight it and its counterpart
  obj = get(gcf,'currentobj');
  ho = get(findobj(gcf,'tag','orig_sides'),'user');
  hr = get(findobj(gcf,'tag','rect_sides'),'user');
  k = find((obj==ho) | (obj==hr));
  if ~isempty(k)
    set(ho(k),'linesty',':')
    set(hr(k),'linesty',':')
    drawnow
    % When button is released, we'll revert to previous state
    set(gcf,'windowbuttonup','cback(''sideup'');')
  end
  
case 'sideup'
  % End side highlighting and restore vertex selection mode
  ho = get(findobj(gcf,'tag','orig_sides'),'user');
  hr = get(findobj(gcf,'tag','rect_sides'),'user');
  set(ho,'linesty','-')
  set(hr,'linesty','-')
  % Nonobvious trick: Put the vertices back on "top" so they are
  % selectable 
  ax = findobj(gcf,'type','axes');
  set(get(ax(1),'user'),'linesty','.')
  set(get(ax(2),'user'),'linesty','.')
  set(gcf,'windowbuttonup','cback(''select'');')
  
case 'compute'
  % Compute rectified polygon
  betar = get(findobj(gcf,'tag','betar'),'user');
  cr = get(findobj(gcf,'tag','cr'),'user');
  Q = get(findobj(gcf,'tag','Q'),'user');
  % Must be a valid set of angles
  if abs(sum(betar)+2) < 1e-10
    % Do it
    wr = NaN*betar;
    % Reference side will include a nontrivial vertex
    k = min(find(betar<0));
    wr([k rem(k,length(betar(:)))+1]) = [0;1];
    [affr,wr] = craffine(wr,betar,cr,Q);
    set(findobj(gcf,'tag','wr'),'user',wr)
    set(findobj(gcf,'tag','affr'),'user',affr)
    % Update picture
    axes(findobj(gcf,'tag','crrect_rectified'));    
    cla
    axis auto
    [hr,h] = plotpoly(wr,betar,1);
    set(findobj('tag','rect_sides','user',hr(:)))	  
    set(h(:,1),'markersize',20)
    set(h(abs(betar+.5)<eps,:),'color',[1 0 1])
    set(h(abs(betar-.5)<eps,:),'color',[0 .5 .5])
    set(h(abs(betar-1) <eps,:),'color',[0 0 0])
    set(gca,'user',h)
  else
    fprintf('The sum of the beta''s must be -2 to compute the rectified polygon.\n')
  end
  
case 'radioexclude'
  t = findobj(gcf,'tag','angradio');
  k = find( t==gcbo );
  set(t(k),'value',1)
  set(t([1:k-1 k+1:4]),'value',0)
end

