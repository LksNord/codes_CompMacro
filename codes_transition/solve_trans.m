function [TRANS] = solve_trans(param,grid,tau,A,SS_init,SS_final)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to solve the transition of an Aiyagari model
% HH problem solved by EGM
% inputs:
%       - param: structure containing the necessary parameter values
%       - grid: structure containing grids
%               (z needs to be specified, k will be created)
%       - tau: path for labor income tax rate
%       - SS_init: initial steady state
%       - SS_final: final steady state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Preallocation
    r                   = [linspace(SS_init.r,SS_final.r,param.TT/2), SS_final.r*ones(1,param.TT/2)] ; % !!! this is the initial guess for the path of r !!!
    dist                = zeros(param.nz,param.nkap,param.TT);
    cpol                = zeros(param.nz,param.nkap,param.TT);
    kpol                = zeros(param.nz,param.nkap,param.TT);
    Ksup                = zeros(1,param.TT);
    C                   = zeros(1,param.TT);

% initialize final and initial conditions (taken as given, will not be updated across iterations!)
    dist(:,:,1)         = SS_init.dist;
    kpol(:,:,end)       = SS_final.kpol;
    cpol(:,:,end)       = SS_final.cpol;
    Ksup(1)             = sum(sum(dist(:,:,1).*(ones(param.nz,1)*grid.k)));

% iterate on transition path
    for iter=1:param.maxiter
        
        % given the guess of interest rate recover equilibrium objects
            w = (1-param.alpha)*(A.*(param.alpha./(r+param.delta)).^param.alpha).^(1/(1-param.alpha));
            T = tau.*w.*param.labor;

        % Solve policy functions (iterate backwards)
            for t=param.TT-1:-1:1

                % solve one time step of the HH problem
                    [cpol(:,:,t),kpol(:,:,t)] = stepEGM(cpol(:,:,t+1),r(t),r(t+1),w(t),tau(t),T(t),param,grid);
         
            end


        % Compute distribution (iterate forward)
        
            dist(:,:,2:end) = 0; % reset to zeros (important as cumulative method!)
            C(1)                = sum(sum(dist(:,:,1).*cpol(:,:,1))); % first period consumption

            for t=1:param.TT-1

                % iterate distribution forward
                    [dist(:,:,t+1)] = stepDist(dist(:,:,t),kpol(:,:,t),param,grid);
                    
                % compute HH aggregates
                    Ksup(t+1) = sum(sum(dist(:,:,t).*kpol(:,:,t)));
                    C(t+1) = sum(sum(dist(:,:,t+1).*cpol(:,:,t+1)));
                     
            end


        % implied interest rate
            r_new           = A.*param.alpha.*Ksup.^(param.alpha-1)*param.labor^(1-param.alpha)-param.delta;

        % display progress
            disp( [ 'iter = ', num2str(iter),'   ', 'err_r = ', num2str((max(abs(r - r_new))))])

        % Check if new interest rate is the same as old one, otherwise go
        % upgrade interest rate guess and do the previous steps again
        if (max(abs(r - r_new)) < param.tol_r)
            break
        else
            if (iter<=10) % begin with slower updating for stability
                weight = param.weightr/2;
            else
                weight = param.weightr;
            end
            r = (1-weight)*r + weight*r_new;
        end


    end

    if iter == param.maxiter
        disp(['Maximum number of iterations on transition reached. Remaining error is ', num2str(max(abs(r - r_new)))])
    end
% load output
    TRANS.cpol = cpol;
    TRANS.kpol = kpol;
    TRANS.dist = dist;
    TRANS.r = r;
    TRANS.tau = tau;
    TRANS.T = T;
    TRANS.w = w;
    TRANS.Ksup = Ksup;
    TRANS.C = C;
    TRANS.Y = A.*Ksup.^param.alpha.*param.labor^(1-param.alpha);




