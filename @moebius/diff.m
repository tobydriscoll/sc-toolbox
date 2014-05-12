function f = diff(map)
%Differentiate a Moebius transformation.
%   DIFF(M) returns a callable function that evaluates to the derivative of 
%   the Moebius transformation M.

%   Copyright (c) 2007 by Toby Driscoll.

a = map.coeff(1); b = map.coeff(2); c = map.coeff(3); d = map.coeff(4);
f = @(z) (b*c-a*d)./(d*z + c).^2;