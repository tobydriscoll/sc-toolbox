function acc = accuracy(M)
%ACCURACY Apparent accuracy of Schwarz-Christoffel CR rectified map.

%   Copyright 1998 by Toby Driscoll.
%   $Id: accuracy.m 7 1998-05-10 04:37:19Z tad $

% Just return accuracy of underlying crdiskmap.
acc = accuracy(M.diskmap);
