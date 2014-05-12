function wq_sum = wqsum(wa,phia,kwa,ic,wb,phib,radius,u,w0,w1,nptq,qwork,linearc,data,param4)
%   WQSUM calculates the complex integral from wa to wb along a line
%   segment or a circular arc with possible singularity at wa. The 
%   definition of WQSUM can be found in USER'S GUIDE to DSCPACK, Chenglie Hu 1995

%   Index arrangement:
    iwt1=nptq*(ic*data.M+kwa-1)+1;
    if (kwa == 0) 
        iwt1=nptq*(data.M+data.N)+1;
    end
    iwt2=iwt1+nptq-1;
    ioffst=nptq*(data.M+data.N+1);

%   Compute Gauss-Jacobi sum(w(j)*prod(x(j))):
    if (linearc==1)
%   Integrate along a circular arc:
        pwh = (phib-phia)/2;
        pwc = (phib+phia)/2;
        uw0 = u*w0;
        u_w1 = u ./ w1;
        w = radius*exp(i*(pwc+pwh*qwork(iwt1:iwt2)));
        wq_sum = i*pwh*qwork(ioffst+iwt1:ioffst+iwt2).*w*wprod(w,u,uw0,u_w1,data,param4);
        return;
    else
%   Integrate along a line segment:
      wh = (wb-wa)/2;
      wc = (wa+wb)/2;
      uw0 = u*w0;
      u_w1 = u ./ w1;
      w = wc+wh*qwork(iwt1:iwt2);
      wq_sum = wh*qwork(ioffst+iwt1:ioffst+iwt2)*wprod(w,u,uw0,u_w1,data,param4);     
      return;
    end