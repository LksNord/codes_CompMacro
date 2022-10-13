function [J,Jsup,Jdem] = compute_Jac(param,grid,SS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to combe the sequence space Jacobian of an AIyagari model for
% a productivity shock. HH problem solved by EGM
% inputs:
%       - param: structure containing the necessary parameter values
%       - grid: structure containing grids
%               (z needs to be specified, k will be created)
%       - SS: steady state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 0 - Some Initialization
    dr = 0.00000001; % size of deviation to which we are computing a numerical derivative

%% Step 1 - derivative of policies
% make one backward iteration to shock in T-1 to get derivative of HH policies
    
    % interest rate path with minor deviation
        r = SS.r*ones(1,param.TT);
        r(param.TT-1)=SS.r+dr;
        
    % Preallocation
        cpol                = zeros(param.nz,param.nkap,param.TT);
        kpol                = zeros(param.nz,param.nkap,param.TT);

    % initialize final and initial conditions (taken as given, will not be updated across iterations!)
        kpol(:,:,end)       = SS.kpol;
        cpol(:,:,end)       = SS.cpol;


    % given the interest rate path recover the wage
        %w = SS.w*ones(1,param.TT);
        w = (1-param.alpha)*(param.A.*(param.alpha./(r+param.delta)).^param.alpha).^(1/(1-param.alpha));

    % Solve policy functions (iterate backwards)
        for t=param.TT-1:-1:1

            % solve one time step of the HH problem
                [cpol(:,:,t),kpol(:,:,t)] = stepEGM(cpol(:,:,t+1),r(t),r(t+1),w(t),0,0,param,grid);

        end

%% Step 2 - Compute time zero derivatives dK_0s and dD_1s

    % initialize
        dK_0s = NaN(1,param.TT-1);
        dD_1s = NaN(param.nz*param.nkap,param.TT-1);
        [~,LambdaSS] = getDist_continuous(param,grid,SS.kpol,1);
    % loop over all possible dates of the shock s
        for t = param.TT-1:-1:1
            
            s = param.TT-t; % distance from today to the shock

            % compute dK_0s
            dK_0s(s) = (sum(sum(SS.dist.*kpol(:,:,t)))-sum(sum(SS.dist.*SS.kpol)))/dr;

            % compute dD_1s
            [~,Lambda0] = getDist_continuous(param,grid,kpol(:,:,t),1);
            dD_1s(:,s) = (reshape(SS.dist,1,param.nz*param.nkap)*(Lambda0-LambdaSS))/dr;

        end

%% Step 3 - Compute expectation vectors E_t's
        
    % initialize
        E_t=NaN(param.nz*param.nkap,param.TT-2);
        
    % fill for all t
        Lambda_rec = (LambdaSS^0);
        y_SS = reshape(SS.kpol,param.nz*param.nkap,1);
        for t = 1:param.TT-2            
            E_t(:,t) = Lambda_rec * y_SS;
            Lambda_rec = Lambda_rec*LambdaSS;
        end

%% Step 4 - Create F matrix
        
    % initialize
        F = NaN(param.TT-1, param.TT-1);

    % first row
        F(1,:)=dK_0s;

    % other rows
        for s=1:param.TT-1
            for t=2:param.TT-1
                F(t,s)= E_t(:,t-1)'*dD_1s(:,s);
            end
        end
 
 %% Step 5 - Create Jacobian of Capital Supply
        
    % initialize
        Jsup = NaN(param.TT-1, param.TT-1);

    % first row / column
        Jsup(1,:)=F(1,:);
        Jsup(:,1)=F(:,1);

    % other rows
        for s=2:param.TT-1
            for t=2:param.TT-1
                Jsup(t,s)= Jsup(t-1,s-1)+F(t,s);
            end
        end
    
 %% Step 6 - Create Jacobian of Capital Demand
    
    % derivative of Kdem at the steady state interest rate
        dKdem = (1/(param.alpha-1)) * (1/(param.alpha*param.A*param.labor^(1-param.alpha))) * ((SS.r+param.delta)/(param.alpha*param.A*param.labor^(1-param.alpha)))^((2-param.alpha)/(param.alpha-1));
    
    % Jacobian
        Jdem = dKdem * eye(param.TT-1);
 
 %% Step 7 - Create Jacobian of Capital Market Clearing
 
    J = [zeros(1,param.TT-1); Jsup(1:end-1,:)] - Jdem;
    
    % careful here about the timing of the supply Jacobian. row t is
    % savings in t, i.e. capital supply in t+1