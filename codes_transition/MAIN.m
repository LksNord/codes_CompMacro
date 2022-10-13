%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solving Transitions of the Aiyagari model                 %%
%% collected for Alex Monge's Quant Macro Course @ EUI       %%
%% by Lukas Nord, November 2021                              %%
%%                                                           %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These codes solve for transitions between two steady states of a baseline
% Aiyagari model. The code can produce either a deterministic transition
% to a new steady state after a change in the labor tax rate or the IRF of
% the economy in response to a one-off (MIT) productivity shock.

% The codes are meant to illustrate solution methods as easily
% accessible as possible and are not necessarily optimized for maximum
% speed. Suggestions to further improve the codes are very welcome.

% Please direct errors or questions to lukas.nord@eui.eu.

clear;
close all;
clc;


%% Initialisation

% algorithm
    param.tol_r           = 0.000001;          % tolerance for convergence of interest rate
    param.tol_pol         = 0.000001;         % tolerance for convergence of policy functions
    param.TT              = 200;             % length of transition path
    param.maxiter         = 100;             % maximum number of iterations on transition path
    param.weightr         = 0.1;             % updating weight for r transition path

% parameters
    param.gamma           = 2;               % risk aversion 
    param.beta            = 0.94;            % subjective discount factor 
    param.delta           = 0.1;            % depreciation
    param.A               = 1;               % aggregate productivity
    param.alpha           = 0.36;            % capital's share of income
    param.nkap            = 200;             % number of asset grid points
    param.b               = -2;              % exogenous borrowing limit
    param.B               = 20;             % upper bounnd on assets
    grid.k                = linspace(param.b,param.B,param.nkap); % equally spaced grid
    %grid.k                = param.b + linspace(0,1,param.nkap).^2*(param.B-param.b); % grid more dense at lower end, increase exponents to increase mass at lower end

% for the chebyshev approximation
    param.Ncheb = 12;                        % order of the polynomials considered
    param.Mcheb = 20;                       % number of collocated points
    
% discretizing AR(1) for income   
    % the process we approximate: log(s_t) = rho*log(s_t-1)+sigmaint*sqrt(1-rho^2)*error_t
    param.nz              = 5;                % number of discretized income states
    param.rho             = 0.9;              % first-order autoregressive coefficient of income
    param.sigmaLR         = 0.2;              % long-run standard deviation of the income process
    param.sigma           = param.sigmaLR*sqrt(1-param.rho^2); % standard deviation of error_t
    
    %Pz is transition matrix of the Markov chain
    %logz is the discretized states of log labor earnings
    %distz is the invariant distribution of Markov chain
    [grid.Pz,grid.logz,grid.distz] = markovappr(param.rho,param.sigma,3,param.nz);
    grid.z               = exp(grid.logz);    % bring back to levels
    param.labor          = grid.z*grid.distz; % aggregate labor is average efficiency units   

%% Shocks
% tax system
    tau_init = 0.1;                     % labor income tax in initial steady state
    tau_final = 0.1;%0.25;%                 % labor income tax in final steady state
    tau_path = [linspace(tau_init,tau_final,35), tau_final*ones(1,param.TT-35)]; % transition path for the tax rate

% productivity
    A_path = param.A*ones(1,param.TT);
    A_path(1) = 0.95*param.A;%param.A;%
    for t = 2:50 % careful to chose long enough for convergence to original level
        A_path(t) = (1-0.8)*param.A+0.8*A_path(t-1); % mean reversion in productivity
    end
        
%% Solve for the initial Steady State
    [SS_init] = solve_aiyagari(param,grid,tau_init);
    fprintf('Initial steady state computed. \n\n')

%% Solve for the final Steady State
    [SS_final] = solve_aiyagari(param,grid,tau_final);
    fprintf('Final steady state computed. \n\n')

%% Solve transition
    [TRANS] = solve_trans(param,grid,tau_path,A_path,SS_init,SS_final);
    fprintf('Transitional dynamics solved. \n\n')

%% Plot Paths

    figure(1)
    plot(-4:1:param.TT,[ones(1,5)*SS_init.r,TRANS.r],'LineWidth',2)
    xlim([-4,param.TT])
    xlabel('t','FontSize',14)
    ylabel('r','FontSize',14)
    title('interest rate path','FontSize',14)

    figure(2)
    plot(-4:1:param.TT,[ones(1,5)*SS_init.k,TRANS.Ksup],'LineWidth',2)  
    xlim([-4,param.TT])
    xlabel('t','FontSize',14)
    ylabel('K','FontSize',14)
    title('capital path','FontSize',14)
    
    figure(3)
    plot(-4:1:param.TT,[ones(1,5)*SS_init.T,TRANS.T],-4:1:param.TT,[ones(1,5)*SS_init.tau,TRANS.tau],'LineWidth',2)  
    xlim([-4,param.TT])
    xlabel('t','FontSize',14)
    title('tax system path','FontSize',14)
    legend({'T','\tau'},'FontSize',12,'Location','SouthEast')

    figure(4)
    plot(-4:1:param.TT,[ones(1,5)*SS_init.w,TRANS.w],'LineWidth',2)
    xlim([-4,param.TT])
    xlabel('t','FontSize',14)
    ylabel('w','FontSize',14)
    title('wage path','FontSize',14)