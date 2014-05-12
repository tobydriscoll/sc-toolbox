function w_dsc = wdsc(zz,u,c,w0,w1,phi0,phi1,nptq,qwork,eps,iopt,dataz)
%   WDSC computes the inverse map after all necessary parameters have
%   been determined. Euler's scheme for solving ODE is used to give the
%   initial guess which is then refined by newton's iteration.
%   zz is the point at which inverse map is evaluated. eps is the required
%   accuracy supplied by the user.
%   Note: zz is not allowed to be a vertex!

    param4 = thdata(u);
    uw0 = u*w0; % needed by wprod
    u_w1 = u ./ w1;% needed by wprod
    
%  Generate the initial value for solving ODE.
   [knz,inz] = nearz(zz,dataz);
  
   if (inz==2)
       disp('The Inverse Evaluation Failed')
       exit;
   end
   if (inz == 0)
        if (knz >= 2)
            wi = (w0(knz)+w0(knz-1))/2;
        else
            wi = (w0(knz)+w0(dataz.M))/2;
        end
    else
        if (knz >= 2)
            wi = (w1(knz)+w1(knz-1))/2;
        else 
            wi = (w0(knz)+w1(dataz.N))/2;
        end
            wi = u*wi/abs(wi);
    end 

    zs = zdsc(wi,0,2,u,c,w0,w1,phi0,phi1,nptq,qwork,1,dataz);

%  Solve ODE initial value problem (along line segment from zs to zz)
%  by Euler's scheme to generate the initial guess:
        
    for K = 1:20
          wi = wi + (zz-zs)/ (20*c*wprod(wi,u,uw0,u_w1,dataz,param4));
    end
    if (abs(wi) > 1)
        wi = wi/ (abs(wi)+abs(1-abs(wi)));
    end
    if (abs(wi) < u) 
        wi = u*wi*(abs(wi)+abs(u-abs(wi)))/abs(wi);
    end

%  Refine the solution by Newton's iteration:
    IT = 1;
    wfn = zz - zdsc(wi,0,2,u,c,w0,w1,phi0,phi1,nptq,qwork,iopt,dataz);
    while((abs(wfn) >= eps) && (IT <=15) )
          wi = wi + wfn/(c*wprod(wi,u,uw0,u_w1,dataz,param4));
          if (abs(wi) > 1) 
              wi = wi/(abs(wi)+abs(1-abs(wi)));
          end
          if (abs(wi) < u)
              wi = u*wi* (abs(wi)+abs(u-abs(wi)))/abs(wi);
          end
          IT = IT + 1;
          wfn = zz - zdsc(wi,0,2,u,c,w0,w1,phi0,phi1,nptq,qwork,iopt,dataz);
    end

    if (abs(wfn) < eps)
        w_dsc = wi;
        % Verification:
        zz1 = zdsc(wi,0,2,u,c,w0,w1,phi0,phi1,nptq,qwork,1,dataz);  
    else
        %  The iteration failed to meet the tolerance in 15 iterations
        %  Try a different vertex as a reference point:
        zz1 = complex(1,1) + zz;
    end

    if ((abs(wi) >= u) && (abs(wi) <= 1) && (abs(zz-zz1) <= 1e-3))
        return;
    end
    
    temp = dataz;
    if (inz == 0)
          temp.Z0(knz) = zz;
    else
          temp.Z1(knz) = zz;
    end 
    
    %Recursive steps
    w_dsc = wdsc(zz,u,c,w0,w1,phi0,phi1,nptq,qwork,eps,iopt,temp);
 
    return;