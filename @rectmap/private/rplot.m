function [H,RE,IM] = rplot(w,beta,z,c,L,re,im,options)
%RPLOT  Image of cartesian grid under Schwarz-Christoffel rectangle map.
%   RPLOT(W,BETA,Z,C,L) will adaptively plot the images under the
%   Schwarz-Christoffel rectangle map of ten evenly spaced horizontal
%   and vertical lines in the retangle RECT.  The arguments are as in
%   RPARAM.
%
%   RPLOT(W,BETA,Z,C,L,M,N) will plot images of M evenly spaced vertical
%   and N evenly spaced horizontal lines.
%
%   RPLOT(W,BETA,Z,C,L,RE,IM) will plot images of vertical lines whose
%   real parts are given in RE and horizontal lines whose imaginary
%   parts are given in IM.  Either argument may be empty.
%
%   RPLOT(W,BETA,Z,C,L,RE,IM,OPTIONS) allows customization of RPLOT's
%   behavior.  See SCPLTOPT.
%
%   H = RPLOT(W,BETA,Z,C,L,...) returns a vector of handles to all the
%   curves drawn in the interior of the polygon.  [H,RE,IM] =
%   RPLOT(W,BETA,Z,C,L,...) also returns the abscissae and ordinates of
%   the lines comprising the grid.
%
%   See also SCPLTOPT, RPARAM, RMAP, RDISP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: rplot.m 298 2009-09-15 14:36:37Z driscoll $

w = w(:);
beta = beta(:);
z = z(:);
[w,beta,z,corners] = rcorners(w,beta,z);
rect = z(corners);
import sctool.*

% Parse input
if nargin < 8
    options = [];
    if nargin < 7
        im = [];
        if nargin < 7
            re = [];
        end
    end
end

Kp = imag(rect(2));
K = rect(1);

% Empty arguments default to 10
if isempty([re(:);im(:)])
    re = 10;
    im = 10;
end

% Integer arguments must be converted to specific values
if (length(re)==1) && (re == round(re))
    if re < 1
        re = [];
    else
        m = re;
        re = linspace(-K,K,m+2);
        re([1,m+2]) = [];
    end
end
if (length(im)==1) && (im == round(im))
    if im < 1
        im = [];
    else
        m = im;
        im = linspace(0,Kp,m+2);
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

% Plot vertical lines...
linh = gobjects(length(re),2);
for j = 1:length(re)
    % Start evenly spaced
    zp = re(j) + 1i*linspace(0,Kp,15).';
    new = true(size(zp));
    wp = NaN(length(zp),1);
    
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
    
    % Adaptive refinement to make smooth curve
    iter = 0;
    while (any(new)) && (iter < maxrefn)
        drawnow
        neww = rmap(zp(new),w,beta,z,c,L,qdat);
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
        [zp,wp,new] = scpadapt(zp,wp,minlen,maxlen,axlim);
    end
    % Set the lines to be solid
    if verLessThan('matlab','8.4')
        set(linh(j,1),'erasemode','back')
        set(linh(j,1),'marker','none','linestyle','-','user',zp)
        if draw2
            % Replace the points with the endpoints
            set(linh(j,2),'erasemode','back')
            set(linh(j,2),'marker','none','linestyle','-',...
                'xdata',re(j)*[1 1],'ydata',[0 Kp])
        end
    else
        clearpoints(linh(j,1))
        addpoints(linh(j,1),real(wp(~isnan(wp))),imag(wp(~isnan(wp))));
        set(linh(j,1),'marker','none','linestyle','-','user',zp)
        if draw2
            % Replace the points with smooth curve
            clearpoints(linh(j,2))
            addpoints(linh(j,2),re(j)*[1 1],[0 Kp])
            set(linh(j,2),'marker','none','linestyle','-')
        end       
    end
    drawnow
end

% Plot horizontal lines...
linh1 = linh;
linh = gobjects(length(im),2);
for j = 1:length(im)
    % Start evenly spaced
    zp = linspace(-K,K,15).' + 1i*im(j);
    new = true(size(zp));
    wp = NaN(length(zp),1);
    
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
    
    % Adaptive refinement to make smooth curve
    iter = 0;
    while (any(new)) && (iter < maxrefn)
        drawnow
        neww = rmap(zp(new),w,beta,z,c,L,qdat);
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
        [zp,wp,new] = scpadapt(zp,wp,minlen,maxlen,axlim);
    end
    % Set the lines to be solid
    if verLessThan('matlab','8.4')
        set(linh(j,1),'erasemode','back')
        set(linh(j,1),'marker','none','linestyle','-','user',zp)
        if draw2
            % Replace the points with the endpoints
            set(linh(j,2),'erasemode','back')
            set(linh(j,2),'marker','none','linestyle','-',...
                'xdata',[-K K],'ydata',im(j)*[1 1])
        end
    else
        clearpoints(linh(j,1))
        addpoints(linh(j,1),real(wp(~isnan(wp))),imag(wp(~isnan(wp))));
        set(linh(j,1),'marker','none','linestyle','-','user',zp)
        if draw2
            % Replace the points with (hopefully) a smooth circle
            clearpoints(linh(j,2))
            addpoints(linh(j,2),[-K K],im(j)*[1 1])
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
        RE = re;
        if nargout > 2
            IM = im;
        end
    end
end

end
