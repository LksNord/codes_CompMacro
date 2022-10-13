%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solving Transitions of the Aiyagari model                 %%
%% collected for Alex Monge's Quant Macro Course @ EUI       %%
%% by Lukas Nord, November 2021                              %%
%%                                                           %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These codes solve for the dynamics of an Aiyagari model in response to a 
% one-off (MIT) productivity shock. The household problem is solved by
% endogenous grid method. For the aggregate dynamics we apply:
% a) extended path algortihm (as previously)
% b) updates of the interest path using the sequence space Jacobian
%       ( see Auclert et al. (ECMA,2021))

% The codes are meant to illustrate solution methods as easily
% accessible as possible and are not necessarily optimized for maximum
% speed. Suggestions to further improve the codes are very welcome.

% Please direct errors or questions to lukas.nord@eui.eu.

clear;
close all;
clc;


%% Initialisation

% algorithm
    param.tol_r           = 1e-6;          % tolerance for convergence of interest rate
    param.tol_pol         = 1e-12;         % tolerance for convergence of policy functions
    param.TT              = 250;             % length of transition path
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

%% Productivity Shock
    A_path = param.A*ones(1,param.TT);
    A_path(1) = 0.95*param.A;%param.A;%
    for t = 2:50 % careful to chose long enough for convergence to original level
        A_path(t) = (1-0.8)*param.A+0.8*A_path(t-1); % mean reversion in productivity
    end
        
%% Solve for the Steady State
    [SS] = solve_aiyagari(param,grid,0);
    fprintf('Initial steady state computed. \n\n')

%% Solve transition with extended path
    tic
    [TRANS] = solve_trans(param,grid,A_path,SS,SS);
    TRANS.time = toc;
    fprintf('Transitional dynamics solved, it took me %3.4f seconds. \n\n',[TRANS.time])

%% Solve transition with sequence space Jacobian
    
    tic;
    J = compute_Jac(param,grid,SS); % J is the Jacobian of the capital market clearing condition
    timeJ = toc;
    fprintf('Jacobian computed, it took me %3.4f seconds. \n\n',[timeJ])

    
    tic
    [TRANS_SSJ] = solve_transSSJ(param,grid,A_path,SS,J);
    TRANS_SSJ.time = toc;
    fprintf('Transitional dynamics (SSJ) solved, it took me %3.4f seconds. \n\n',[TRANS_SSJ.time])


%% Plot Paths
    
    % differences in paths between methods
        maxdiff_r = max(abs(TRANS.r-TRANS_SSJ.r));
        maxdiff_K = max(abs(TRANS.Ksup-TRANS_SSJ.Ksup));
        maxdiff_w = max(abs(TRANS.w-TRANS_SSJ.w));

    
    figure(1)
    plot(-4:1:param.TT,[ones(1,5)*SS.r,TRANS.r],-4:1:param.TT,[ones(1,5)*SS.r,TRANS_SSJ.r],'LineWidth',2)
    xlim([-4,param.TT])
    xlabel('t','FontSize',14)
    ylabel('r','FontSize',14)
    title('interest rate path','FontSize',14)
    legend({'trans','trans_{SSJ}'},'FontSize',12,'Location','SouthEast')


    figure(2)
    plot(-4:1:param.TT,[ones(1,5)*SS.k,TRANS.Ksup],-4:1:param.TT,[ones(1,5)*SS.k,TRANS_SSJ.Ksup],'LineWidth',2)
    xlim([-4,param.TT])
    xlabel('t','FontSize',14)
    ylabel('K','FontSize',14)
    title('capital path','FontSize',14)
    legend({'trans','trans_{SSJ}'},'FontSize',12,'Location','SouthEast')

    
    figure(3)
    plot(-4:1:param.TT,[ones(1,5)*SS.w,TRANS.w],-4:1:param.TT,[ones(1,5)*SS.w,TRANS_SSJ.w],'LineWidth',2)
    xlim([-4,param.TT])
    xlabel('t','FontSize',14)
    ylabel('w','FontSize',14)
    title('wage path','FontSize',14)
    legend({'trans','trans_{SSJ}'},'FontSize',12,'Location','SouthEast')
