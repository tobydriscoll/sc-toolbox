function w = dscinvmap(z,map)
% DSCINVMAP Map a pt z to a pt w

% DSCINVMAP(z,map) maps a pt z in the doubly connected region to the
% annulus. map is an annulusmap object.

% see @annulusmap/annulusmap.m

% check if z is in the doubly connected region
% Not implemented yet.
%----------------------
%----------------------

% z is not allowed to be a vertex.
% check if z is in Z0.
idx = min(find(map.Z0 == z));
if isempty(idx)==0
    error('The point calculated is a vertex.');
else
    % check if z is in Z1.
    idx = min(find(map.Z1 == z));
    if isempty(idx)==0
       error('The point calculated is a vertex.');
    end
end
nptq = 8;
eps = 1e-9;

%Making the bridge to old subroutine
dataz = struct('M',map.M,'N',map.N,'Z0',map.Z0,'Z1',map.Z1,'ALFA0',map.ALFA0,'ALFA1',map.ALFA1,'ISHAPE',map.ISHAPE);
w = wdsc(z,map.u,map.c,map.w0,map.w1,map.phi0,map.phi1,nptq,map.qwork,eps,1,dataz);