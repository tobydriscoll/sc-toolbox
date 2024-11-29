function qwork = qinit(dataz,nptq)
% QINIT computes the Gauss-Jacobi nodes & weights for Gauss-Jacobi quadrature.
% Work array qwork must be dimensioned 2*nptq*(M+N+1). It is divided into 2*(M+N+1) 
% vectors of length nptq: The first M+N+1 contain quadrature nodes on output, the next
% M+N+1 contain the corresponding weights on output. nptq is the number of
% G-J nodes (same as weights) used.

import sctool.gaussj;

%   For each finite vertex, compute nodes & weights
%   For one-sided Gauss-Jacobi quadrature 

import sctool.gaussj
for K=1:(dataz.M+dataz.N)
    inodes = nptq*(K-1)+1;
    iwts = nptq*(dataz.M+dataz.N+K)+1;
    if (K <= dataz.M)
        alpha = dataz.ALFA0(K)-1;
        if (dataz.ALFA0(K) > 0)
         [qwork(inodes:inodes+nptq-1),qwork(iwts:iwts+nptq-1)]=gaussj(nptq,0,alpha);                  
        else
         qwork(iwts:iwts+nptq-1)=0;
         qwork(inodes:inodes+nptq-1)=0;              
        end
    else
        alpha = dataz.ALFA1(K-dataz.M)-1;
        [qwork(inodes:inodes+nptq-1),qwork(iwts:iwts+nptq-1)]=gaussj(nptq,0,alpha);
    end

% Take singularities into account in advance for the purpose of saving
% certain amount of calculation in WQSUM
    qwork(iwts:iwts+nptq-1)=qwork(iwts:iwts+nptq-1).*((1+qwork(inodes:inodes+nptq-1)).^(-alpha));
end
% Compute nodes & weights for pure Gaussian quadrature:
    inodes = nptq*(dataz.M+dataz.N)+1;
    iwts = nptq*(2*(dataz.M+dataz.N)+1)+1;
    [qwork(inodes:inodes+nptq-1),qwork(iwts:iwts+nptq-1)]=gaussj(nptq,0,0);
