function [u,c,w0,w1,phi0,phi1] = dscsolv(iguess,nptq,qwork,ishape,linearc,dataz)
%  DSCSOLV solves the nonlinear system for D-SC parameters.
%  The fortran subroutine of DSCSOLV can be found in USER'S GUIDE to 
%  DSCPACK, Chenglie Hu 1995

% Provide space if ishape = 0
w1 = zeros(1,dataz.N);
phi1 = zeros(1,dataz.N);
nshape = 0;
%ind(1) = 0;
ind = []; %modified
if (ishape==1)
% Determine the number of outer boundary components,etc:
    nshape = 1;
        for I = 2:dataz.M - 1
            if (dataz.ALFA0(I) <= 0)
                nshape = nshape + 1;
                ind(nshape) = I;
            end
        end
    ind(1) = 0;    
    nshape = nshape + 1;
    ind(nshape) = dataz.M - 1;
end
  
% fix one prevertex
w0(dataz.M) = 1;
phi0(dataz.M) = 0;
x(2)=0;
x(3)=0;
% INITIAL GUESS (iguess=0):
if (iguess==0)
    x(1) = 1/0.5 - 1/0.46;
    ave = 2*pi/dataz.N;
    x(4:dataz.N+1) = log((ave + 0.0001*(1:dataz.N-2))./(ave + 0.0001*(2:dataz.N-1)));
    x(dataz.N + 2) = log((ave + 0.0001*(dataz.N-1))/(2*pi-(dataz.N-1)*(ave + (dataz.N)*0.00005)));
    x(dataz.N + 3) = 1/(4-0.1) - 1/(4+0.1);
%    x(dataz.N + 3) = 1/(pi-0.1) - 1/(pi+0.1);
%    x(dataz.N + 3) = 1/(pi-0.001) - 1/(pi+0.001);
    ave = 2*pi/(dataz.M);
    x(dataz.N+4:dataz.M+dataz.N+1) = log((ave+0.0001*(0:dataz.M-3))./(ave + 0.0001*(1:dataz.M-2)));
    x(dataz.M + dataz.N + 2) = log((ave+0.0001*(dataz.M-2))/(2*pi-(dataz.M-1)*(ave+(dataz.M-2)*0.00005)));
elseif (iguess==1)
% INITIAL GUESS (iguess=1):
    % Initial guess can be modified for faster convergence
    x(1) = 1/0.53 - 1/0.43; 
    x(4:dataz.N+2) = 0;
    x(dataz.N + 3) = 1/(4-0.1) - 1/(4+0.1);
%    x(dataz.N + 3) = 1/(pi-0.001) - 1/(pi+0.001);
    x(dataz.N+4:dataz.M+dataz.N+2) = 0;
else 
    x(1) = 1/(0.98-u) - 1/(u-0.02);
    for K = 1:dataz.N - 1
        KN = K - 1;
        if (KN == 0) 
            KN = dataz.N;
        end
            top = phi1(K) - phi1(KN);
        if (top < 0) 
            top = top + 2*pi;
        end
            botm = phi1(K+1) - phi1(K);
        if (botm < 0) 
            botm = botm + 2*pi;
        end
            x(3+K) = log(top) - log(botm);
    end 
        x(dataz.N + 3) = 1/(pi-phi1(dataz.N)) - 1/(pi+phi1(dataz.N));
        for K = 1:dataz.M - 1
            KM = K - 1;
            if (KM==0) 
                KM = dataz.M;
            end
                top = phi0(K) - phi0(KM);
            if (top < 0) 
                top = top + 2*pi;
            end
                botm = phi0(K+1) - phi0(K);    
            if (botm < 0) 
                botm = botm + 2*pi;
            end
              x(dataz.N + 3 + K) = log(top) - log(botm);
          end 
end

% Calculate the initial guess x(2) & x(3) to match
% the choice of x(1),x(4),...,x(M+N+2):
[u,c,w0,w1,phi0,phi1] = xwtran(x,w0,w1,phi0,phi1,dataz);
param4 = thdata(u);
wint = wquad(w0(dataz.M),0,dataz.M,0,w1(dataz.N),0,dataz.N,1,0,u,w0,w1,nptq,qwork,0,2,dataz,param4);
c1 = (dataz.Z1(dataz.N)-dataz.Z0(dataz.M))/wint;
x(2) = real(c1);
x(3) = imag(c1);

% nesolve,fsolve,SDogleg passing variables:
fdat = {w0,w1,phi0,phi1,nptq,qwork,linearc,ishape,nshape,ind,dataz};

% Solve nonlinear systems by using SDogleg:
%disp('Wait..be patient :) ')
%[x,info] = sdogleg(@dscfun,fdat,x,[1 1e-9 1e-7 1e-7 1000 1e-3 ]);

% Solve nonlinear systems by using nesolve:
%tol = 1e-8;
%opt = zeros(16,1);
%opt(1) = 2;
%opt(2) = 2;
%opt(6) = 100*(16-3);
%opt(8) = tol;
%opt(9) = min(eps^(2/3),tol/10);
%opt(12) = nptq;
%x = nesolve(@dscfun,x(:),opt,fdat);
%x = nesolve(@dscfun,x(:),2,fdat); %old style

% Solve nonlinear systems by using fsolve:
options=optimset('Display','iter');   % Option to display output
x = fsolve(@dscfun,x',options,fdat);

% Convert back x to DSC parameters u,c,w0,w1,phi0,phi1
[u,c,w0,w1,phi0,phi1] = xwtran(x,w0,w1,phi0,phi1,dataz);
return;
