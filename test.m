clear;
IPOLY = 5;

if IPOLY==1
    q = sqrt(2);
    a = 1+q;
    z0 = [a+a*i -a+a*i -a-a*i a-a*i];
    z1 = [q q*i -q -q*i];
    Q = sqrt(2);
    p0 = polygon(z0);
    p1 = polygon(z1);
    rgn = dscpolygons(p1,p0)
elseif IPOLY==2
    % Use sdogleg instead. nesolve is very slow
    z0 = [0.5+2.5*i 0.5+0.5*i 1+0.5*i 1+i 1+0.5*i 0.5+0.5*i 0.5+2.5*i 2.5*i 0 2 2+2.5*i];
    z1 = [1+2*i 1+1.5*i 1+2*i 1.5+2*i 1.5+0.5*i 1.5+2*i];
    p0 = polygon(z0);
    p1 = polygon(z1);
    rgn = dscpolygons(p1,p0)
elseif IPOLY==3
    z0 = [-1-i 1-i 1+i -1+i];
    z1 = [0.5*i 0 -0.5 0 -0.5*i -0.5-0.5*i 0.5-0.5*i 0.25-0.5*i ... 
          0.25-0.25*i 0.25-0.5*i -0.5*i 0 0.5 0 0.5*i 0.5+0.5*i -0.5+0.5*i];
    p0 = polygon(z0);
    p1 = polygon(z1);
    rgn = dscpolygons(p1,p0)
elseif IPOLY==4
    z0 = [-2-i 2-i 2+2*i -0.8+2*i 1+0.5*i -1+2*i -2+2*i];
    z1 = [0 -1];
    p0 = polygon(z0);
    p1 = polygon(z1);
    rgn = dscpolygons(p1,p0)
elseif IPOLY == 5
    z0 = [1+i 4+i 4+4*i 1+4*i];
    z1 = [2.75+1.5*i 3.25+3.25*i 1.5+3*i];
    p0 = polygon(z0);
    p1 = polygon(z1);
    rgn = dscpolygons(p1,p0)
end

tic
f = annulusmap(rgn)
toc
tic
plot(f);
toc