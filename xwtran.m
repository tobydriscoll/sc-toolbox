function [u,c,w0,w1,phi0,phi1] = xwtran(x,w0,w1,phi0,phi1,dataz)
%  XWTRAN transforms x(k) (unconstrained parameters suggested by Daeppen)
%  to actual DSC parameters : u,c,w0,w1. phi0 & phi1 are arguments of the
%  prevertices contained in w0 & w1.

if (abs(x(1)) <= 1e-14)
    u = 0.50;
else
    u = (x(1)-2-sqrt(0.9216*x(1)^2+4))/(2*x(1));
    u = (0.0196*x(1)-1)/(u*x(1));
end

c = complex(x(2),x(3));
if (abs(x(dataz.N + 3)) <= 1e-14)
    phi1(dataz.N) = 0;
else
    ph = (1+sqrt(1+(pi^2)*(x(dataz.N + 3)^2)))/x(dataz.N + 3);
    phi1(dataz.N) = (pi^2)/ph;
end
dph = 1;
phsum = dph;
for k = 1:dataz.N - 1
    dph = dph/exp(x(3+k));
    phsum = phsum + dph;
end
dph = 2*pi/phsum;
phi1(1) = phi1(dataz.N) + dph;
w1(1) = u*complex(cos(phi1(1)),sin(phi1(1)));
w1(dataz.N) = u*complex(cos(phi1(dataz.N)),sin(phi1(dataz.N)));
phsum = phi1(1);
for k = 1:dataz.N - 2
    dph = dph/exp(x(3+k));
    phsum = phsum + dph;
    phi1(k+1) = phsum;
    w1(k+1) = u*complex(cos(phsum),sin(phsum));
end
dph = 1;
phsum = dph;
for k = 1:dataz.M - 1
    dph = dph/exp(x(dataz.N + 3 + k));
    phsum = phsum + dph;
end
dph = 2*pi/phsum;
phsum = dph;
phi0(1) = dph;
w0(1) = complex(cos(dph),sin(dph));
for k = 1:dataz.M - 2
    dph = dph/exp(x(dataz.N + 3 + k));
    phsum = phsum + dph;
    phi0(k+1) = phsum;
    w0(k+1) = complex(cos(phsum),sin(phsum));
end
