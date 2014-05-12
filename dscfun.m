function fval = dscfun(x,fdat)
% DSCFUN forms the nonlinear system satisfied by d-sc parameters.
% This function will be called by nonlinear system solver nesolve.
% The fortran subroutine of DSCFUN can be found in USER'S GUIDE to 
% DSCPACK, Chenglie Hu 1995
[w0,w1,phi0,phi1,nptq,qwork,linearc,ishape,nshape,ind,dataz] = deal(fdat{:});
%fval = zeros(dataz.M + dataz.N + 2,1);

% Transform unconstrained x(k) to actual d-sc parameters
% and compute dataz to be used in wtheta: 
[u,c,w0,w1,phi0,phi1] = xwtran(x,w0,w1,phi0,phi1,dataz);
param4 = thdata(u);

% Two equations to eliminate possible rotation of the inner polygon:
win1 = dataz.Z1(1)-dataz.Z1(dataz.N)-c*wquad(w1(dataz.N),phi1(dataz.N),dataz.N,1,w1(1),phi1(1),1,1,u,u,w0,w1,nptq,qwork,linearc,2,dataz,param4);
fval(1) = real(win1);
fval(2) = imag(win1);

%  N-1 side length conditions for the inner polygon
for I = 1:dataz.N - 1
wint1 = wquad(w1(I),phi1(I),I,1,w1(I+1),phi1(I+1),I+1,1,u,u,w0,w1,nptq,qwork,linearc,2,dataz,param4);
fval(I+2) = abs(dataz.Z1(I+1)-dataz.Z1(I)) - abs(c*wint1);
end

% Two equations to fix the relative position of the inner polygon:
test1 = cos(phi1(dataz.N));

if (test1 >= u)
    win2 = wquad(w0(dataz.M),0,dataz.M,0,w1(dataz.N),0,dataz.N,1,0,u,w0,w1,nptq,qwork,0,2,dataz,param4);    
else

% If the line path from w0(M) to w1(N) is out of domain, the combination
% of two different paths will be used instead:
    wx = u;
    wline = wquad(w0(dataz.M),0,dataz.M,0,wx,0,0,2,0,u,w0,w1,nptq,qwork,0,2,dataz,param4);
      if (phi1(dataz.N) <= 0)
          warc = wquad(w1(dataz.N),phi1(dataz.N),dataz.N,1,wx,0,0,2,u,u,w0,w1,nptq,qwork,1,2,dataz,param4);
          win2 = wline - warc;
      else
          warc = wquad(wx,0,0,2,w1(dataz.N),phi1(dataz.N),dataz.N,1,u,u,w0,w1,nptq,qwork,1,2,dataz,param4);
          win2 = wline + warc;
      end    
end

fval(dataz.N + 2) = real(dataz.Z1(dataz.N)-dataz.Z0(dataz.M)-c*win2);
fval(dataz.N + 3) = imag(dataz.Z1(dataz.N)-dataz.Z0(dataz.M)-c*win2);

% Two equations to eliminate possible rotation of the outer polygon:
win3 = dataz.Z0(1)-dataz.Z0(dataz.M)-c*wquad(w0(dataz.M),phi0(dataz.M),dataz.M,0,w0(1),phi0(1),1,0,1,u,w0,w1,nptq,qwork,linearc,2,dataz,param4);
fval(dataz.N + 4) = real(win3);
fval(dataz.N + 5) = imag(win3);

if (dataz.M == 3)
    fval = fval(:); %MMMMM
    return;
end

if (ishape==1)
    ;
else
% M-3 side length conditions of the outer polygon:
for J = 1:dataz.M - 3
    wint2 = wquad(w0(J),phi0(J),J,0,w0(J+1),phi0(J+1),J+1,0,1,u,w0,w1,nptq,qwork,linearc,2,dataz,param4);
    fval(dataz.N + 5 + J) = abs(dataz.Z0(J+1)-dataz.Z0(J)) - abs(c*wint2);
end  
fval = fval(:); %MMMMM
return;
end

% Outer polygon contains some infinite vertices & for each of them
% two length conditions will be replaced by a complex integral:
for K = 1:nshape - 1
    if ((ind(K+1)==2) | (ind(K) >= (ind(K+1)-2)))
        ;
    else
        for J = (ind(K) + 1):(ind(K+1) - 2)
          wint3 = wquad(w0(J),phi0(J),J,0,w0(J+1),phi0(J+1),J+1,0,1,u,w0,w1,nptq,qwork,linearc,2,dataz,param4);
          fval(dataz.N + 5 + J) = abs(dataz.Z0(J+1)-dataz.Z0(J)) - abs(c*wint3);
        end 
    end
  
   if ((K==(nshape-1)) | (ind(K+1) == (dataz.M - 1)))
       ;
   else
    % The combination fo three different paths is used to integrate from wa
    % to wb to avoid domain program. The by-product of this is that it is
    % numerically more efficient than a single line path:
          wa = w0(ind(K+1)-1);
          wb = w0(ind(K+1)+1);
          phia = phi0(ind(K+1)-1);
          phib = phi0(ind(K+1)+1);
          radius = (1+u)/2;
          wai = radius*exp(i*phia);
          wbi = radius*exp(i*phib);
          wline1 = wquad(wa,0,ind(K+1)-1,0,wai,0,0,2,0,u,w0,w1,nptq,qwork,0,2,dataz,param4);
          wline2 = wquad(wb,0,ind(K+1)+1,0,wbi,0,0,2,0,u,w0,w1,nptq,qwork,0,2,dataz,param4);
          wcircle = wquad(wai,phia,0,2,wbi,phib,0,2,radius,u,w0,w1,nptq,qwork,1,2,dataz,param4);
          win4 = c*(wline1+wcircle-wline2);
          fval(dataz.N+5+ind(K+1)-1) = real(dataz.Z0(ind(K+1)+1)-dataz.Z0(ind(K+1)-1)-win4);
          fval(dataz.N+5+ind(K+1)) = imag(dataz.Z0(ind(K+1)+1)-dataz.Z0(ind(K+1)-1)-win4);  
   end
end
fval = fval(:); %MMMMM
return;