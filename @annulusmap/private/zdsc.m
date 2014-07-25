function z_dsc = zdsc(ww,kww,ic,u,c,w0,w1,phi0,phi1,nptq,qwork,iopt,dataz)
% ZDSC computes the forward map z(ww) by integrating from wa to ww, where
% wa is the nearest prevertex to ww. kww=k, ic = 0 if ww = w0(k), or kww =
% k, ic=1 if ww=w1(k), or kww = 0 and ic=2 otherwise.
% iopt: = 1 is normally assumed for any geometry.
%       = # other than 1 assumes that only line segment path is used.

param4 = thdata(u);
ibd = 1;
ww0 = ww;

if (iopt ~= 1) 
    ;
else
    %if (abs(abs(ww0)-1) <= 1e-11) %original
    if (abs(abs(ww0)-1) == 0) 
        ww0 = (1+u)*ww0/2;
        ibd = 2;
    end
end

% Find the nearest prevertex to the pt ww0
[knear,inear] = nearw(ww0,w0,w1,dataz);
if (inear==0)
    za = dataz.Z0(knear);
    wa = w0(knear);
else
    za = dataz.Z1(knear);
    wa = w1(knear);
end

if (iopt~=1)
    z_dsc = za + c*wquad(wa,0,knear,inear,ww,0,kww,ic,0,u,w0,w1,nptq,qwork,0,1,dataz,param4);
    return;
else

% Determine the point closest to wa on the same concentric circle as ww0 is
% on:
    if (imag(ww0)>=0)
        phiww0 = angle(ww0);
    else
        phiww0 = angle(ww0)+2*pi;
    end

    dww0 = abs(ww0);
    if (inear==0)
        phiwb = phi0(knear);
    else
        phiwb = phi1(knear);
    end 
    wb = dww0*exp(i*phiwb);

% Integration from wa to wb on a line segment:
    %if (abs(wb-wa) <= 1e-11) %original
    if (abs(wb-wa) == 0)
        wint1 = 0;
    else
         wint1 = wquad(wa,0,knear,inear,wb,0,0,2,0,u,w0,w1,nptq,qwork,0,1,dataz,param4); 
    end

% Integration from wb to ww on a circular arc:
    %if (abs(wb-ww0) <= 1e-11) %original
    if (abs(wb-ww0) == 0)
          wint2 = 0;
          % Evaluate the mapping function. The integration path is
          % a combination of a circular arc and line segment(s):
          z_dsc = za + c* (wint1+wint2);
          if (ibd==2)
          z_dsc = z_dsc + c*wquad(ww0,0,0,2,ww,0,kww,ic,0,u,w0,w1,nptq,qwork,0,1,dataz,param4);
          end
      return;
    end

    if (abs(phiwb-2*pi-phiww0) < abs(phiwb-phiww0)) 
        phiwb = phiwb - 2*pi;
    end 
    if (abs(phiww0-2*pi-phiwb) < abs(phiwb-phiww0)) 
        phiww0 = phiww0 - 2*pi;
    end
      
    %if (abs(wb-wa) <= 1e-11) %original
    if (abs(wb-wa) == 0)   
        wint2 = wquad(wb,phiwb,knear,inear,ww0,phiww0,kww,ic,dww0,u,w0,w1,nptq,qwork,1,1,dataz,param4);
        % Evaluate the mapping function. The integration path is
        % a combination of a circular arc and line segment(s):
        z_dsc = za + c*(wint1+wint2);
        if (ibd==2)
          z_dsc = z_dsc + c*wquad(ww0,0,0,2,ww,0,kww,ic,0,u,w0,w1,nptq,qwork,0,1,dataz,param4);
        end
      return;
    end

    wint2 = wquad(wb,phiwb,0,2,ww0,phiww0,0,2,dww0,u,w0,w1,nptq,qwork,1,1,dataz,param4);

 % Evaluate the mapping function. The integration path is
 % a combination of a circular arc and line segment(s):   
    z_dsc = za + c*(wint1+wint2);
    if (ibd==2)
        z_dsc = z_dsc + c*wquad(ww0,0,0,2,ww,0,kww,ic,0,u,w0,w1,nptq,qwork,0,1,dataz,param4);
    end
    return; 
end