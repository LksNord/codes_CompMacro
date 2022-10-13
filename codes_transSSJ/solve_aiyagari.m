function [SS] = solve_aiyagari(param,grid,tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to solve the GE of an Aiyagari model
% the HH problem is solved by EGM
% inputs:
%       - param: structure containing the necessary parameter values
%       - grid: structure containing grids
%               (z needs to be specified, k will be created)
%       - tau: labor income tax rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialisation
    minrate     = -param.delta;
    maxrate     = (1-param.beta)/(param.beta);  
    err         = 1;
    iter        = 0;
    
% Main Aiyagari loop to iterate on equilibrium interest rate
    while abs(err)>param.tol_r
        iter = iter + 1;
        
    % update guesses
        r0      = 0.5*(maxrate+minrate);                                   % update interest rate guess (bisection)
        k0      = ((r0+param.delta)/(param.alpha*param.A*param.labor^(1-param.alpha)))^(1/(param.alpha-1));    % implied guess for capital stock
        w0      = (1-param.alpha)*(param.A*(param.alpha/(r0+param.delta))^param.alpha)^(1/(1-param.alpha));    % wage implied by interest rate guess
        T       = tau*w0*param.labor;                                      % implied transfer to clear government budget
        
    % test whether natural borrowing limit holds
        if (r0>0) && (param.b < -((1-tau)*w0*grid.z(1)+T)/r0)
            fprintf('Warning! Natural borrowing limit violated!\n')
        end

    
     % solve HH problem for given interest rate
        [cpol,kpol] = solveHH_EGM(r0,w0,tau,T,param,grid);
     
     % compute distribution over households
        [dist] = getDist_continuous(param,grid,kpol,0);
        
     % implied aggregate capital supplied by households
        k1 = sum(sum(dist.*kpol));
        r1 = param.alpha*param.A*param.labor^(1-param.alpha)*max(k1,0.01)^(param.alpha-1)-param.delta;
        
     % update interest rate guess
        err = r1-r0;
        if err < 0
            maxrate = r0;
        else
            minrate = r0;
        end
        disp(['k0 = ',num2str(k0),', k1 = ',num2str(k1),', r0 = ',num2str(r0),', r1 = ',num2str(r1)])

    end
    
% load output structure
    SS.r = r0;
    SS.w = w0;
    SS.k = k0;
    SS.cpol = cpol;
    SS.kpol = kpol;
    SS.dist = dist;
    SS.gridk = grid.k;
    SS.tau = tau;
    SS.T = T;
    SS.y = param.A*param.labor^(1-param.alpha)*k1^param.alpha;