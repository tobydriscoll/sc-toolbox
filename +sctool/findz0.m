function [z0,w0] = findz0(prefix,wp,mapfun,w,beta,z,c,qdat,aux)
%SCIMAPZ0 (not intended for calling directly by the user)
%   SCIMAPZ0 returns starting points for computing inverses of
%   Schwarz-Christoffel maps.
%
%   Each wp(j) (in the polygon plane) requires z0(j) (in the fundamental
%   domain) whose image w0(j) is such that the line segment from w0(j)
%   to wp(j) lies in the target (interior or exterior) region.  The
%   algorithm here is to choose z0(j) as a (weighted) average of
%   successive pairs of adjacent prevertices.  The resulting w0(j) is on
%   a polygon side.  Each choice is tested by looking for intersections
%   of the segment with (other) sides of the polygon.
%
%   After randomly trying 10 weights with such prevertex pairs, the
%   routine gives up.  Failures are pretty rare.  Slits are the most
%   likely cause of trouble, since the intersection method doesn't know
%   "which side" of the slit it's on.  In such a case you will have to
%   supply starting points manually, perhaps by a continuation method.
%
%   See also HPINVMAP, DINVMAP, DEINVMAP, RINVMAP, STINVMAP.

%   Copyright 1997 by Toby Driscoll.  Last updated 05/07/97.

%   P.S. This file illustrates why the different domains in the SC
%   Toolbox have mostly independent M-files.  The contingencies for the
%   various geometries become rather cumbersome.

n = length(w);
shape = wp;
wp = wp(:);
z0 = wp;
w0 = wp;
from_disk = strcmp(prefix(1),'d');
from_hp = strcmp(prefix,'hp');
from_strip = strcmp(prefix,'st');
from_rect = strcmp(prefix,'r');

% Calculate arguments of the directed polygon edges.
if from_strip
    % Renumber to put left end of strip first
    atinf = find(isinf(z));
    renum = [atinf(1):n 1:atinf(1)-1];
    w = w(renum);
    z = z(renum);
    beta = beta(renum);
    qdat(:,1:n) = qdat(:,renum);
    qdat(:,n+1+(1:n)) = qdat(:,n+1+renum);
    kinf = find(isinf(z), 1, 'last' );
    argw = cumsum([angle(w(3)-w(2));-pi*beta([3:n,1])]);
    argw = argw([n,1:n-1]);
else
    argw = cumsum([angle(w(2)-w(1));-pi*beta(2:n)]);
end

% Express each side in a form to allow detection of intersections.
infty = isinf(w);
fwd = [2:n 1];
anchor = zeros(1,n);
anchor(~infty) = w(~infty);
anchor(infty) = w( fwd(infty) );        % use the finite endpoint
direcn = exp(1i*argw);
direcn(infty) = -direcn(infty);         % reverse
len = abs( w(fwd) - w );

if from_disk
    argz = angle(z);
    argz(argz<=0) = argz(argz<=0) + 2*pi;
end

if from_rect
    % Extra argument given
    L = qdat;
    qdat = aux;
end

factor = 0.5;				% images of midpoints of preverts
done = zeros(1,length(wp));
m = length(wp);
iter = 0;
if length(qdat) > 1
    tol = 1000*10^(-size(qdat,1));
else
    tol = qdat;
end

zbase = NaN*ones(n,1);
wbase = NaN*ones(n,1);
idx = [];
while m > 0				% while some not done
    % Choose a "basis point" on each side of the polygon.
    for j = 1:n
        if from_disk
            if j<n
                zbase(j) = exp(1i*(factor*argz(j) + (1-factor)*argz(j+1)));
            else
                zbase(j) = exp(1i*(factor*argz(n) + (1-factor)*(2*pi+argz(1))));
            end
        elseif from_hp
            if j < n-1			% between two finite points
                zbase(j) = z(j) + factor*(z(j+1)-z(j));
            elseif j==n-1			% between x(n-1) & Inf
                zbase(j) = max(10,z(n-1))/factor;
            else				% between -Inf and x(1)
                zbase(j) = min(-10,z(1))/factor;
            end
        elseif from_strip
            if j==1
                zbase(j) = min(-1,real(z(2)))/factor;
            elseif j==kinf-1
                zbase(j) = max(1,real(z(kinf-1)))/factor;
            elseif j==kinf
                zbase(j) = 1i+max(1,real(z(kinf+1)))/factor;
            elseif j==n
                zbase(j) = 1i+min(-1,real(z(n)))/factor;
            else
                zbase(j) = z(j) + factor*(z(j+1)-z(j));
            end
        elseif from_rect
            zbase(j) = z(j) + factor*(z(rem(j,n)+1)-z(j));
            % Can't use 0 or iK' as basis points.
            if abs(zbase(j)) < 1e-4
                zbase(j) = zbase(j) + .2i;
            elseif abs(zbase(j)-1i*max(imag(z))) < 1e-4
                zbase(j) = zbase(j) - .2i;
            end
        end
        
        % Find images of all the z-plane basis points.
        wbase(j) = feval(mapfun,zbase(j));
        
        % Project each base point exactly to the nearest polygon side. This
        % is needed to make intersection detection go smoothly in borderline
        % cases.
        proj = real( (wbase(j)-anchor(j)) * conj(direcn(j)) );
        wbase(j) = anchor(j) + proj*direcn(j);
        
    end
    
    if isempty(idx)
        % First time through, assign nearest basis point to each image point
        [dist,idx] = min(abs( wp(~done,ones(n,1)).' - wbase(:,ones(m,1)) ));
    else
        % Other times, just change those that failed.
        idx(~done) = rem(idx(~done),n) + 1;
    end
    z0(~done) = zbase(idx(~done));
    w0(~done) = wbase(idx(~done));
    
    % Now, cycle thru basis points
    for j = 1:n
        % Those points who come from basis j and need checking
        active = (idx==j) & (~done);
        if any(active)
            % Assume for now that it's good
            done(active) = ones(1,sum(active));
            % Test line segment for intersections with other sides.
            % We'll parameterize line segment and polygon side, compute parameters
            % at intersection, and check parameters at intersection.
            for k=[1:j-1,j+1:n]
                A(:,1) = [ real(direcn(k)); imag(direcn(k)) ];
                for p = find(active)
                    dif = (w0(p)-wp(p));
                    A(:,2) = [real(dif);imag(dif)];
                    % Get line segment and side parameters at intersection.
                    if rcond(A) < eps
                        % All 4 involved points are collinear.
                        wpx = real( (wp(p)-anchor(k)) / direcn(k) );
                        w0x = real( (w0(p)-anchor(k)) / direcn(k) );
                        if (wpx*w0x < 0) || ((wpx-len(k))*(w0x-len(k)) < 0)
                            % Intersection interior to segment: it's no good
                            done(p) = 0;
                        end
                    else
                        dif = (w0(p)-anchor(k));
                        s = A\[real(dif);imag(dif)];
                        % Intersection occurs interior to side? and segment?
                        if s(1)>=0 && s(1)<=len(k)
                            if abs(s(2)-1) < tol
                                % Special case: wp(p) is on polygon side k
                                z0(p) = zbase(k);
                                w0(p) = wbase(k);
                            elseif abs(s(2)) < tol
                                % Happens when two sides are partially coincident (slit)
                                % Check against normal to that side
                                if real( conj(wp(p)-w0(p))*1i*direcn(k) ) > 0
                                    % Line segment came from "outside"
                                    done(p) = 0;
                                end
                            elseif s(2) > 0 && s(2) < 1
                                % Intersection interior to segment: it's no good
                                done(p) = 0;
                            end
                        end
                    end
                end
            end
            
            % Short circuit if everything is finished
            m = sum(~done);
            if ~m, break, end
        end
    end
    if iter > 2*n
        error('Can''t seem to choose starting points.  Supply them manually.')
    else
        iter = iter + 1;
    end
    factor = rand(1);			% abandon midpoints
end

shape(:) = z0;
z0 = shape;
shape(:) = w0;
w0 = shape;
