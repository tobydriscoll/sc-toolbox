function w = evalinv(map,z)
%EVALINV   Evaluate the inverse map.
%   EVALINV(MAP,Z) finds the inverse image of a point Z under the annulusmap 
%   MAP. That is, it maps from a doubly connected polygonal domain to the
%   canonical annulus.
%
% See also ANNULUSMAP, ANNULUSMAP.EvAL.

% Copyright by Toby Driscoll, 2014.
% Written by Alfa Heryudono, 2003. 

% TODO: check if z is in the doubly connected region

% z is not allowed to be a vertex.
% check if z is in Z0.
idx = find(map.Z0 == z, 1 );
if isempty(idx)==0
    error('The point calculated is a vertex.');
else
    % check if z is in Z1.
    idx = find(map.Z1 == z, 1 );
    if isempty(idx)==0
       error('The point calculated is a vertex.');
    end
end
nptq = 8;
eps = 1e-9;

%Making the bridge to old subroutine
dataz = struct('M',map.M,'N',map.N,'Z0',map.Z0,'Z1',map.Z1,'ALFA0',map.ALFA0,'ALFA1',map.ALFA1,'ISHAPE',map.ISHAPE);
w = wdsc(z,map.u,map.c,map.w0,map.w1,map.phi0,map.phi1,nptq,map.qwork,eps,1,dataz);