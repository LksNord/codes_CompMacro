function [cheby_vec] = cheby_vec(x,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constructs vector of chebyshev polynomials up to order N evaluated at x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    cheby_vec=NaN(N+1,1);
    
    cheby_vec(1) = 1;
    cheby_vec(2) = x;
    
    for i=2:(N)
        cheby_vec(i+1) = 2*x*cheby_vec(i) - cheby_vec(i-1);
    end
    
    

    
