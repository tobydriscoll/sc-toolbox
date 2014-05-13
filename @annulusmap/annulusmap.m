function map = annulusmap(rgn)
%ANNULUSMAP Solve the nonlinear system for D-SC parameter.
%   ANNULUSMAP(RGN) solves the nonlinear system for D-SC parameter of a
%   doubly connected region RGN as created by DSCPOLYGONS.
%
%   See also DSCPOLYGONS.

%   Copyright by Alfa Heryudono, 2003.

% Get inner and outer polygons from region rgn
[p1 p0] = get(rgn);

% Making the bridge to old subroutine :-) (I am lazy to clean the code
% right now)
map.region = rgn; %store rgn. Needed by plot. Too many unnecessary assignment. clean code later.

% Check if the outer polygon is unbounded (ISHAPE = 1)
if isinf(p0)==1
    % Truncate Infinite vertices.
    p0 = truncate(p0);
    map.M = length(p0); map.N = length(p1);
    map.Z0 = vertex(p0).'; map.Z1 = vertex(Z1).';
    map.ALFA0 = angle(p0)'; map.ALFA1 = 2 - angle(p1)';
    map.ISHAPE = 1;
else
    p0 = truncate(p0);
    map.M = length(p0); map.N = length(p1);
    map.Z0 = vertex(p0).'; map.Z1 = vertex(p1).';
    map.ALFA0 = angle(p0)'; map.ALFA1 = 2 - angle(p1)';
    map.ISHAPE = 0;
end

% Generate the Gauss-Jacobi weights & nodes
nptq = 8;
map.qwork = qinit(map,nptq);

% Solve D-SC parameter
iguess = 0;  %(0 nonequally spaced guess or 1 equally spaced guess)
linearc = 1; % (0 line path or 1 circular path)
[map.u,map.c,map.w0,map.w1,map.phi0,map.phi1] = dscsolv(iguess,nptq,map.qwork,map.ISHAPE,linearc,map);

map = class(map,'annulusmap');