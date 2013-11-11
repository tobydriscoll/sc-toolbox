function [H,RE,IM] = crrplot(w,beta,wr,betar,cr,aff,affr,Q,re,im,options)
%CRRPLOT  Image of cartesian grid under Schwarz-Christoffel rectified map.
%   CRRPLOT(W,BETA,WR,BETAR,CR,AFF,AFFR,Q) will adaptively plot the
%   images under the Schwarz-Christoffel rectified map of ten evenly
%   spaced horizontal and vertical lines in the polygon WR.
%
%   With additional arguments M and N, there will be M evenly spaced
%   vertical and N evenly spaced horizontal lines. If these arguments
%   are instead vectors (or noninteger scalars) RE and IM, the vertical
%   lines will have real parts given by RE and the horizontal lines will
%   have imaginary parts given by IM. Either argument may be empty.
%
%   CRRPLOT(W,BETA,WR,BETAR,CR,AFF,AFFR,Q,M,N,OPTIONS) allows
%   customization of CRRPLOT's behavior.  See SCPLTOPT.
%
%   H = CRRPLOT(...) returns a vector of handles to all the curves drawn
%   in the interior of the polygon W.  [H,RE,IM] = CRRPLOT(...) also
%   returns the abscissae and ordinates of the lines comprising the
%   grid.
%       
%   See also SCPLTOPT, CRPARAM, CRRECT, CRRMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: crrplot.m 226 2003-01-08 16:59:17Z driscoll $

% Parse input and initialize
n = length(w);
w = w(:);
beta = beta(:);
wr = wr(:);
betar = betar(:);

if nargin < 11
  options = [];
  if nargin < 10
    im = [];
    if nargin < 9
      re = [];
    end
  end
end

xlim = [min(real(wr)) max(real(wr))];
ylim = [min(imag(wr)) max(imag(wr))];

% Empty arguments default to 10
if isempty([re(:);im(:)])
  re = 10;
  im = 10;
end

% Integer arguments must be converted to specific values
if (length(re)==1) & (re == round(re))
  if re < 1
    re = [];
  else
    m = re;
    re = linspace(xlim(1),xlim(2),m+2);
    re([1,m+2]) = [];
  end
end
if (length(im)==1) & (im == round(im))
  if im < 1
    im = [];
  else
    m = im;
    im = linspace(ylim(1),ylim(2),m+2);
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
len = max(diff(get(ax(1),'xlim')),diff(get(ax(1),'ylim')));
minlen = len*minlen;
maxlen = len*maxlen;
qdat = scqdata(beta,nqpts);
qdatr = scqdata(betar,nqpts);
tol = 10^(-nqpts);
axlim = axis;

% Do it
hv = crrplot0('ver',re,minlen,maxlen,maxrefn,tol,...
    w,beta,wr,betar,cr,aff,affr,Q,qdat,qdatr,ax,axlim,vis,draw2);
hh = crrplot0('hor',im,minlen,maxlen,maxrefn,tol,...
    w,beta,wr,betar,cr,aff,affr,Q,qdat,qdatr,ax,axlim,vis,draw2);
linh = [hv;hh];
if ~draw2
  linh = linh(:,1);
end
set(linh,'erasemode','normal')

% Force redraw to get clipping enforced
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


function linh = crrplot0(direcn,val,minlen,maxlen,maxrefn,tol,...
                         w,beta,wr,betar,cr,aff,affr,Q,qdat,qdatr, ...
			 ax,axlim,vis,draw2)

%   This is the routine that actually draws the curves corresponding
%   to either vertical and horizontal lines. 

n = length(w);

% Find size of the rectified domain.
siz = max( max(imag(wr))-min(imag(wr)), max(real(wr))-min(real(wr)) );

linh = zeros(length(val),2);
for j = 1:length(val)

  % Find intersections with rectified polygon sides
  if strcmp(direcn(1:3),'ver')
    d = real(wr)-val(j);
    cross = find(diff(sign(d([1:n 1]))) | abs(d) < tol );
    [srtcross,idx] = sort(imag(wr(cross)));
  else    
    d = imag(wr)-val(j);
    cross = find(diff(sign(d([1:n 1]))) | abs(d) < tol );
    [srtcross,idx] = sort(real(wr(cross)));
  end
  cross = cross(idx);

  % Remove (near-) duplicates
  identical = find( abs(diff(srtcross)) < tol );
  srtcross(identical) = [];
  cross(identical) = [];
  numcross = length(cross);

  % Set up intervals between crossings
  interval = [srtcross(1:numcross-1)  srtcross(2:numcross)]';

  % Test points at ends and mindpoints of intervals
  if strcmp(direcn(1:3),'ver')
    testint = val(j) + i*mean(interval);
    testedge = val(j) + i*srtcross;
  else
    testint = i*val(j) + mean(interval);
    testedge = i*val(j) + srtcross;
  end
  
  % Use test points to determine number of curves
  index = isinpoly(testint,wr,betar,tol);
  numcurves = max(isinpoly(testedge,wr,betar,tol));
  numcurves = max(numcurves,max(index));

  % Put initial points in zp (string together intervals).
  delta = diff(interval);
  n0 = max(5, ceil(16*delta/siz));
  zp = interval(1,1);
  for k = 1:length(delta)
    zp = [zp; interval(1,k)+(1:n0(k))'*delta(k)/n0(k)];
  end
  %h = delta/5;
  %zp = interval(ones(6,1),find(index)) + (0:5)'*h(:,find(index));
  %zp = zp(:);
  %zp(find(~diff(zp))) = [];
  
  if strcmp(direcn(1:3),'ver')
    zp = val(j) + i*zp;
  else
    zp = i*val(j) + zp;
  end

  % Prepare for iterative mapping
  lenzp = length(zp);
  new = logical(ones(lenzp,1));
  wp = NaN*ones(lenzp,numcurves);
  qn = NaN*ones(lenzp,numcurves);
  iter = 0;
  
  color = 'k';

  % The individual points will be shown as they are found
  linh(j,1) = line(NaN,NaN,'parent',ax(1),'color',color,'vis',vis{1},...
      'linestyle','none','marker','.','markersize',7,'erasemode','none');
  if draw2
    linh(j,2) = line(NaN,NaN,'parent',ax(2),'color',color,'vis',vis{2},...
	'linestyle','none','marker','.','markersize',7,'erasemode','none');
  end

  % Do the mapping & plotting
  while any(new) & (iter < maxrefn)
    drawnow

    % Map new points 
    [neww,newqn] = crrmap(zp(new),w,beta,wr,betar,cr,aff,affr,...
	Q,qdat,qdatr);
    nc = size(neww,2);
    
    % Check
    if nc > numcurves
      error('Too many values found at some points!')
    end
    
    % Incorporate new points
    wp(new,1:nc) = neww;
    if any(~new)
      qnold = qn;
      qn = NaN*wp;
      qn(~new,:) = qnold;
    end
    qn(new,1:nc) = newqn;

    % Sort the columns so as to make continuous curves
    [wp,qn] = crrsort(wp,qn,Q);

    % Update the points to show progress
    set(linh(j,1),'xdata',real(wp(:)),'ydata',imag(wp(:)))
    if draw2
      set(linh(j,2),'xdata',real(zp),'ydata',imag(zp))
    end

    iter = iter + 1;
    
    % Add points to zp where necessary to make smooth curves
    [zp,wp,new] = scpadapt(zp,wp,minlen,maxlen,axlim);
  end	
  
  % In case we exceeded max # iters, delete the uncomputed points
  wp(new,:) = [];
  
  % Add a row of NaN's. We will stack columns, so the NaN's prevent the
  % graphical joining of separate curves.
  wp(size(wp,1)+1,:) = NaN*wp(1,:);
  
  % Set the lines to be solid
  set(linh(j,1),'erasemode','back')
  set(linh(j,1),'marker','none','linestyle','-')
  set(linh(j,1),'xdata',real(wp(:)),'ydata',imag(wp(:)),'user',zp)
  if draw2
    % Replace the points with the endpoints
    set(linh(j,2),'erasemode','back')
    set(linh(j,2),'marker','none','linestyle','-',...
	'xdata',real(zp([1 end])),'ydata',imag(zp([1 end])))
  end
  drawnow

end


function [wp,qnum] = crrsort(wp,qnum,Q)

%   For multiple-valued maps (self-overlapping source polygons), the
%   CRRMAP function returns columns for all images of each source
%   point. This function sorts those columns for CRRPLOT0 so that
%   each represents a continuous curve.

[len,m] = size(wp);
n3 = size(Q.adjacent,2);

% adj is a matrix of quadrilateral adjacencies. A dummy row and column of
% zeros are added to simulate the "null quadrilateral." The diagonal is set
% to ones to make quads self-adjacent.
adj = zeros(n3+1);
adj(1:n3,1:n3) = Q.adjacent + eye(n3);

% Make the unused positions come from the null quadrilateral
mask = isnan(qnum);
qnum(mask) = (n3+1)*ones(sum(mask(:)),1);

if m > 1
  % For each potential curve
  for k = 1:m
    % Indices of changes in source quadrilateral
    p = find( diff(qnum(1:len,k)) );
    toggle = ~isnan(wp(1,k));
    while ~isempty(p)
      p = p(1);
      % If the change is not to a neighbor, swap columns
      if ~adj(qnum(p,k),qnum(p+1,k))
	% Look for a column that does use a neighbor
	col = find( adj(qnum(p,k),qnum(p+1,:)) );
%%        % Without that, go by physical distance.
%%        if isempty(col)
%%          [tmp,col] = min(abs(wp(p+1,k) - wp(p,:)));
%%          col = col(1);
%%        end
	if ~isempty(col)
	  % Swap
	  wp(p+1:len,[k col]) = wp(p+1:len,[col k]);
	  qnum(p+1:len,[k col]) = qnum(p+1:len,[col k]);
	end
      end
      if toggle
        if isnan(wp(p+1,k)), break, end
      else
        if ~isnan(wp(p+1,k)), toggle=1; end
      end
      % Move to the next change in quads (must recompute)
      p = p + find( diff(qnum(p+1:len,k)) );
    end
  end
end

% Change null quads back to NaNs
mask = (qnum == n3+1);
qnum(mask) = NaN*ones(sum(mask(:)),1);
