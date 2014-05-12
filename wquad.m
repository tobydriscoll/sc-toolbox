function w_quad = wquad(wa,phia,kwa,ica,wb,phib,kwb,icb,radius,u,w0,w1,nptq,qwork,linearc,ievl,dataz,param4)
%   WQUAD calculates the complex integral of WPROD from wa to wb along
%   either a circular arc or a line-segment. Function WQUAD1 is called four
%   times, one for each 1/4 of the interval. Note: WQUAD1 allows only the
%   left endpoint to be a possible singularity. The definition of WQUAD can be 
%   found in USER'S GUIDE to DSCPACK, Chenglie Hu 1995

%   Determine midpts on a line segment or on a circular arc:
if (linearc==0)
    wmid = (wa+wb)/2;
    wmida = (wa+wmid)/2;
    wmidb = (wb+wmid)/2;
    phmid = 0;
    phmida = 0;
    phmidb = 0;
else
    if(ievl ~=1)
        if (phib < phia) 
            phia = phia - 2*pi;
        end
    end
    phmid = (phia+phib)/2;
    wmid = radius*exp(i*phmid);
    phmida = (phia+phmid)/2;
    wmida = radius*exp(i*phmida);
    phmidb = (phib+phmid)/2;
    wmidb = radius*exp(i*phmidb);
end

%   Compound Gauss-Jacobi process according to one-quarter rule:
wa_to_wmida = wquad1(wa,phia,kwa,ica,wmida,phmida,radius,u,w0,w1,nptq,qwork,linearc,dataz,param4);
wmid_to_wmida = wquad1(wmid,phmid,0,2,wmida,phmida,radius,u,w0,w1,nptq,qwork,linearc,dataz,param4);
    wqa = wa_to_wmida - wmid_to_wmida;

wb_to_wmidb = wquad1(wb,phib,kwb,icb,wmidb,phmidb,radius,u,w0,w1,nptq,qwork,linearc,dataz,param4);
wmid_to_wmidb = wquad1(wmid,phmid,0,2,wmidb,phmidb,radius,u,w0,w1,nptq,qwork,linearc,dataz,param4);
    wqb = wb_to_wmidb - wmid_to_wmidb;

    w_quad = wqa - wqb;