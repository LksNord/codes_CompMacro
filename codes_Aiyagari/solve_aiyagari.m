function [SS] = solve_aiyagari(param,grid,method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to solve the GE of an Aiyagari model
% inputs:
%       - param: structure containing the necessary parameter values
%       - grid: structure containing grids
%               (z needs to be specified, k will be created)
%       - method: choice for the solution method applied to the HH problem
%                 1: discrete value function iteration
%                 2: continuous value function iteration
%                 3: policy function iteration
%                 4: endogenous grid method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialisation
    minrate     = -param.delta;
    maxrate     = (1-param.beta)/(param.beta);  
    err         = 1;
    iter        = 0;
    V           = false; % no external guess on first iteration for VFI
    cpol        = false; % no external guess on first iteration for PFI
    
% Main Aiyagari loop to iterate on equilibrium interest rate
    while abs(err)>param.tol_r
        iter = iter + 1;
        
    % update guesses
        r0      = 0.5*(maxrate+minrate);                                   % update interest rate guess (bisection)
        k0      = ((r0+param.delta)/(param.alpha*param.A*param.labor^(1-param.alpha)))^(1/(param.alpha-1));    % implied guess for capital stock
        w0      = (1-param.alpha)*(param.A*(param.alpha/(r0+param.delta))^param.alpha)^(1/(1-param.alpha));    % wage implied by interest rate guess

    % capital grid
        if r0<=0
           blim              = param.b; % imposed borrowing limit if negative interest rate
        else
           blim              = max(param.b, -(w0*grid.z(1))/r0); % check if natural borrowing limit is binding when interest positive          
        end

        grid.k                 = linspace(blim,param.B,param.nkap); % equally spaced grid
        %grid.k                 = blim + linspace(0,1,param.nkap).^2*(param.B-blim); % grid more dense at lower end, increase exponents to increase mass at lower end
    
     % solve HH problem for given interest rate
     if method == 1 % discrete VFI
        [cpol,kpol,kopt,V] = solveHH_VFId(r0,w0,param,grid,V);
     elseif method == 2 % continuous VFI
        [cpol,kpol,V] = solveHH_VFIc(r0,w0,param,grid,V);
     elseif method == 3 % PFI
        [cpol,kpol] = solveHH_PFI(r0,w0,param,grid,cpol);
     elseif method == 4 % EGM
        [cpol,kpol] = solveHH_EGM(r0,w0,param,grid);
     elseif method == 5 % projection  
        [cpol,kpol,V] = solveHH_VFIcheb(r0,w0,param,grid);
     end
     
     % compute distribution over households
     if method == 1
        [dist] = getDist_discrete(param,grid,kopt);
     else
        [dist] = getDist_continuous(param,grid,kpol);
     end
        
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
    
% no value function for Euler Equation based methods
% (can still be computed later from optimal policies)
    if (method==3) || (method==4)
        V = NaN;
    end
    
% load output structure
    SS.r = r0;
    SS.w = w0;
    SS.k = k0;
    SS.cpol = cpol;
    SS.kpol = kpol;
    SS.V = V;
    SS.dist = dist;
    SS.gridk = grid.k;