function w_quad1 = wquad1(wa,phia,kwa,ic,wb,phib,radius,u,w0,w1,nptq,qwork,linearc,dataz,param4)
%   WQUAD1 calculates the complex integral of WPROD from wa to wb along
%   either a circular arc or a line-segment. Compound one-sided
%   Gauss-Jacobi quadrature is used. See subroutine WQSUM for the calling
%   sequence. The definition of WQUAD1 can be found in USER'S GUIDE to DSCPACK, Chenglie Hu 1995

if (abs(wa-wb)==0)
%   Check for the zero length integral
      w_quad1 = 0;
      return;
else
    if (linearc==1) % Circular arc
        % Find minimum distance from wa to the nearest vertex.
        d = [abs(wa-w0),abs(wa-w1)];
        %---------
        % consider all points with distance less than 1e-14 from the vertex
        % are indeed the vertex.
        idx = (d < 1e-15);
        d(idx) = 0;
        %---------
        w_dist = min(2,min(d(find(d))));
        r = min(1,w_dist/abs(wb-wa));
        % -----------------------------
        phaa = phia + r*(phib-phia);
        waa = radius*exp(i*phaa);
        w_quad1 = wqsum(wa,phia,kwa,ic,waa,phaa,radius,u,w0,w1,nptq,qwork,linearc,dataz,param4);
      while (r~=1)
        % Find minimum distance from waa to the nearest vertex.
        d = [abs(waa-w0),abs(waa-w1)];
        %---------
        % consider all points with distance less than 1e-14 from the vertex
        % are indeed the vertex.
        idx = (d < 1e-15);
        d(idx) = 0;
        %---------
        w_dist = min(2,min(d(find(d))));
        r = min(1,w_dist/abs(waa-wb));
        % -----------------------------
        phbb = phaa + r*(phib-phaa);
        wbb = radius*exp(i*phbb);
        w_quad1 = w_quad1 + wqsum(waa,phaa,0,2,wbb,phbb,radius,u,w0,w1,nptq,qwork,linearc,dataz,param4);
        phaa = phbb;
        waa = wbb;
      end
      return;
    else %Line segment
        % Find minimum distance from wa to the nearest vertex.
        d = [abs(wa-w0),abs(wa-w1)];
        %---------
        % consider all points with distance less than 1e-14 from the vertex
        % are indeed the vertex.
        idx = (d < 1e-15);
        d(idx) = 0;
        %---------
        w_dist = min(2,min(d(find(d))));
        r = min(1,w_dist/abs(wb-wa));
        % ------------------------------
        waa = wa + r*(wb-wa);
        w_quad1 = wqsum(wa,0,kwa,ic,waa,0,0,u,w0,w1,nptq,qwork,linearc,dataz,param4);
      while (r~=1)
        % Find minimum distance from waa to the nearest vertex.
        d = [abs(waa-w0),abs(waa-w1)];
        %---------
        % consider all points with distance less than 1e-14 from the vertex
        % are indeed the vertex.
        idx = (d < 1e-15);
        d(idx) = 0;
        %---------
        w_dist = min(2,min(d(find(d))));
        r = min(1,w_dist/abs(waa-wb));
        % -------------------------------
        wbb = waa + r*(wb-waa);
        w_quad1 = w_quad1 + wqsum(waa,0,0,2,wbb,0,0,u,w0,w1,nptq,qwork,linearc,dataz,param4);
        waa = wbb;
      end
      return; 
    end
end