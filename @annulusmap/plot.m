function plot(map,varargin)
%PLOT Graphical representation of an annulusmap.
%   PLOT(MAP) creates a graphical representation of the annulusmap MAP.
%
%   PLOT(map,'num') or PLOT(map,'lab') also plots dots for the vertices and
%   numeric labels of inner polygon and outer polygon. For infinite vertices,
%   two numeric labels are printed.

%   Modification of dplot.
%   Copyright by Alfa Heryudono, 2003.

plotdouble(map.boundary);

hold on

% Making the bridge to the old code
dataz = struct('M',map.M,'N',map.N,'Z0',map.Z0,'Z1',map.Z1,'ALFA0',map.ALFA0,'ALFA1',map.ALFA1,'ISHAPE',map.ISHAPE);

nptq = 8; % Gauss-Jacobi nodes.
param4 = thdata(map.u);


theta = 10;
if (length(theta)==1) && (theta == round(theta))
    m = theta+1;
    theta = linspace(0,2*pi,m);
    theta(m) = [];
end

% Prepare figure
fig = gcf;
figure(fig);
turn_off_hold = ~ishold;

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
color = 'k';

% Plot circles...
R = map.u+(1-map.u)*(1:10)/(11);
linh = zeros(length(R),2);
for j = 1:length(R)
    % Start with evenly spaced theta
    tp = linspace(0,2*pi,100)';
    new = true(length(tp),1);
    wp = NaN*new;
    
    % The individual points will be shown as they are found
    if verLessThan('matlab','8.4')
        % The individual points will be shown as they are found
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
    
    % Adaptively refine theta to make smooth curve
    iter = 0;
    maxrefn = 2;
    while (any(new)) && (iter < maxrefn)
        drawnow
        zp = R(j)*exp(1i*tp(new));
        for k=1:length(zp)
            [knear,inear] = nearw(zp(k),map.w0,map.w1,dataz);
            if(inear == 0)
                zc = map.Z0(knear);
                wa = map.w0(knear);
            else
                zc = map.Z1(knear);
                wa = map.w1(knear);
            end
            neww(k)= zc + map.c*wquad1(wa,0,knear,inear,zp(k),0,0,map.u,map.w0,map.w1,nptq,map.qwork,0,dataz,param4);
        end
        wp(new) = neww;
        iter = iter + 1;
        
        % Update the points to show progress
        if verLessThan('matlab','8.4')
            set(linh(j,1),'xdata',real(wp),'ydata',imag(wp))
            if draw2
                set(linh(j,2),'xdata',R(j)*cos(tp),'ydata',R(j)*sin(tp))
            end
        else
            addpoints(linh(j,1),real(wp(new)),imag(wp(new)))
            if draw2
                addpoints(linh(j,1),R(j)*cos(tp(new)),R(j)*sin(tp(new)))
            end
        end
        drawnow
        
    end
    
    % Set the lines to be solid
    if verLessThan('matlab','8.4')
        set(linh(j,1),'erasemode','back')
        set(linh(j,1),'marker','none','linestyle','-','user',R(j)*exp(1i*tp))
        if draw2
            % Replace the points with (hopefully) a smooth circle
            tp = linspace(0,2*pi,101);
            set(linh(j,2),'erasemode','back')
            set(linh(j,2),'marker','none','linestyle','-',...
                'xdata',R(j)*cos(tp),'ydata',R(j)*sin(tp))
        end
        
    else
        
        clearpoints(linh(j,1))
        addpoints(linh(j,1),real(wp),imag(wp));
        set(linh(j,1),'marker','none','linestyle','-','user',R(j)*exp(1i*tp))
        if draw2
            % Replace the points with (hopefully) a smooth circle
            tp = linspace(0,2*pi,361);
            clearpoints(linh(j,2))
            addpoints(linh(j,2),R(j)*cos(tp),R(j)*sin(tp))
            set(linh(j,2),'marker','none','linestyle','-')
        end
    end
    drawnow
end


neww = [];

% Plot radii...
theta = linspace(0,2*pi,11);
linh1 = linh;
linh = zeros(length(theta),2);
for j = 1:length(theta)
    %Rp = linspace(0,1,14)';
    Rp = linspace(map.u,1,14)';
    zp = Rp*exp(i*theta(j));
    new = true(length(zp),1);
    wp = NaN*new;
    
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
    
    % Adaptively refine to make smooth curve
    iter = 0;
    while (any(new)) & (iter < maxrefn)
        drawnow
        for k=1:length(zp)
            [knear,inear] = nearw(zp(k),map.w0,map.w1,dataz);
            if(inear == 0)
                zc = map.Z0(knear);
                wa = map.w0(knear);
            else
                zc = map.Z1(knear);
                wa = map.w1(knear);
            end
            neww(k)= zc + map.c*wquad1(wa,0,knear,inear,zp(k),0,0,map.u,map.w0,map.w1,nptq,map.qwork,0,dataz,param4);
        end
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
                addpoints(linh(j,1),real(zp(new)),imag(zp(new)))
            end
        end
        
    end
    
    % Set the lines to be solid
    if verLessThan('matlab','8.4')
        set(linh(j,1),'erasemode','back')
        set(linh(j,1),'marker','none','linestyle','-','user',zp)
        if draw2
            % Replace the points with just the ends
            set(linh(j,2),'erasemode','back')
            set(linh(j,2),'marker','none','linestyle','-',...
                'xdata',[0 1]*cos(theta(j)),'ydata',[0 1]*sin(theta(j)))
        end
    else
        clearpoints(linh(j,1))
        addpoints(linh(j,1),real(wp),imag(wp));
        set(linh(j,1),'marker','none','linestyle','-','user',zp)
        if draw2
            % Replace the points with just the ends
            clearpoints(linh(j,2))
            addpoints(linh(j,2),[0 1]*cos(theta(j)),[0 1]*sin(theta(j)))
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

end