function [knear,inear] = nearw(w,w0,w1,dataz)
% NEARW determines the nearest prevertex to a given pt w. On return integer
% inear indicates that the nearest pt is found either in w0 (inear = 0) or
% in w1 (inear = 1). knear contains the corresponding index in w0 or in w1.
% Note: the prevertices corresponding to infinite vertices will be skipped.
% (outer circle only).

dist = 2;
d = abs(w - w0);
knear = find(dataz.ALFA0 > 0 & d < dist);
if isempty(knear)==0
    knear = find(d==min(d(knear)));
    knear = knear(1); % Avoid multiple nearest vertices.
    dist = d(knear);
end

inear = 0;
d = abs(w - w1);
idx = find(d < dist);
if isempty(idx)==0
    knear = find(d==min(d(idx)));
    knear = knear(1); % Avoid multiple nearest vertices.
    inear = 1;
end