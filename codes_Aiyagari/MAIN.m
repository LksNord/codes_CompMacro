%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solving the Aiyagari model with different methods         %%
%% collected for Alex Monge's Quant Macro Course @ EUI       %%
%% by Lukas Nord, November 2021                              %%
%%                                                           %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These codes solve for the steady state of an Aiyagari model,
% using different methods to solve the household problem

% The codes are meant to illustrate different solution methods as easily
% accessible as possible and are not necessarily optimized for maximum
% speed (especially VFIc and PFI). Suggestions to further improve these
% codes are very welcome.

% Please direct errors or questions to lukas.nord@eui.eu.

clear;
close all;
clc;


%% Initialisation

% algorithm
    param.tol_r           = 0.0001;          % tolerance for convergence of interest rate
    param.tol_pol         = 0.00005;         % tolerance for convergence of policy functions
    param.tol_vf          = 0.0001;          % tolerance for convergence of value functions (slower than policies)

% parameters
    param.gamma           = 2;               % risk aversion 
    param.beta            = 0.94;            % subjective discount factor 
    param.delta           = 0.1;            % depreciation
    param.A               = 1;               % aggregate productivity
    param.alpha           = 0.36;            % capital's share of income
    param.nkap            = 200;             % number of asset grid points
    param.b               = -2;              % exogenous borrowing limit
    param.B               = 20;             % upper bounnd on assets

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
    
    
    
%% Solve for the Steady State with discrete VFI
    % solving the Bellman Equation on a discretized asset grid
    tic;
    [SS_VFId] = solve_aiyagari(param,grid,1);
    SS_VFId.time = toc;
    fprintf('Steady state computed (VFId), it took me %3.4f seconds. \n\n',[SS_VFId.time])

%% Solve for the Steady State with continuous VFI
    % solving the Bellman Equation allowing for a continuous asset choice
    % VERY(!) slow, possibly due to incompetent implementation...
    tic;
    [SS_VFIc] = solve_aiyagari(param,grid,2);
    SS_VFIc.time = toc;
    fprintf('Steady state computed (VFId), it took me %3.4f seconds. \n\n',[SS_VFIc.time])

%% Solve for the Steady State with PFI
    % iterating forward on the Euler Equation
    % slow, possibly due to incompetent implementation...
    tic;
    [SS_PFI] = solve_aiyagari(param,grid,3);
    SS_PFI.time = toc;
    fprintf('Steady state computed (PFI), it took me %3.4f seconds. \n\n',[SS_PFI.time])

%% Solve for the Steady State with EGM
    % iterating backward on the Euler Equation
    tic;
    [SS_EGM] = solve_aiyagari(param,grid,4);
    SS_EGM.time = toc;
    fprintf('Steady state computed (EGM), it took me %3.4f seconds. \n\n',[SS_EGM.time])

%% Solve for the Steady State with Projection Method
    % iterating on Bellman Equation approx. with Chebyshev polynomials
    % slow, possibly due to incompetent implementation...
    tic;
    [SS_CHEB] = solve_aiyagari(param,grid,5);
    SS_CHEB.time = toc;
    fprintf('Steady state computed (CHEB), it took me %3.4f seconds. \n\n',[SS_CHEB.time])

%% Plots

    % select for which steady state to plot things
    SS = SS_EGM;

    % Look at some steady state output
    figure(1)
    plot(SS.gridk,SS.cpol(1,:),SS.gridk,SS.cpol(end,:),'LineWidth',2)  
    xlabel('assets','FontSize',14)
    ylabel('consumption','FontSize',14)
    title('consumption policy','FontSize',14)
    legend({'$z_1$','$z_N$'},'FontSize',12,'Location','SouthEast')

    figure(2)
    plot(SS.gridk,SS.kpol(1,:),SS.gridk,SS.kpol(end,:),'LineWidth',2)  
    xlabel('assets','FontSize',14)
    ylabel('savings','FontSize',14)
    title('savings policy','FontSize',14)
    legend({'$z_1$','$z_N$'},'FontSize',12,'Location','SouthEast')

    figure(3)
    yyaxis left
    plot(SS.gridk,sum(SS.dist),'LineWidth',2)
    yyaxis right
    plot(SS.gridk,cumsum(sum(SS.dist)),'LineWidth',2)
    xlabel('assets','FontSize',14)
    title('asset distribution','FontSize',14)
    legend({'PDF (left)','CDF (right)'},'FontSize',12,'Location','EastOutside')


%% Compare methods

    figure(4)
    plot(SS_VFId.gridk,SS_VFId.cpol(1,:),SS_VFIc.gridk,SS_VFIc.cpol(1,:),SS_PFI.gridk,SS_PFI.cpol(1,:),SS_EGM.gridk,SS_EGM.cpol(1,:),SS_CHEB.gridk,SS_CHEB.cpol(1,:),'LineWidth',2)  
    xlabel('assets','FontSize',14)
    ylabel('consumption','FontSize',14)
    title('consumption policy z_1','FontSize',14)
    legend({'VFId','VFIc','PFI','EGM','CHEB'},'FontSize',12,'Location','SouthEast')
    
    figure(5)
    plot(SS_VFId.gridk,SS_VFId.cpol(end,:),SS_VFIc.gridk,SS_VFIc.cpol(end,:),SS_PFI.gridk,SS_PFI.cpol(end,:),SS_EGM.gridk,SS_EGM.cpol(end,:),SS_CHEB.gridk,SS_CHEB.cpol(end,:),'LineWidth',2)  
    xlabel('assets','FontSize',14)
    ylabel('consumption','FontSize',14)
    title('consumption policy z_N','FontSize',14)
    legend({'VFId','VFIc','PFI','EGM','CHEB'},'FontSize',12,'Location','SouthEast')

    figure(6)
    plot(SS_VFId.gridk,sum(SS_VFId.dist),SS_VFIc.gridk,sum(SS_VFIc.dist),SS_PFI.gridk,sum(SS_PFI.dist),SS_EGM.gridk,sum(SS_EGM.dist),SS_CHEB.gridk,sum(SS_CHEB.dist),'LineWidth',2)  
    xlabel('assets','FontSize',14)
    ylabel('PDF','FontSize',14)
    title('asset distribution','FontSize',14)
    legend({'VFId','VFIc','PFI','EGM','CHEB'},'FontSize',12,'Location','SouthEast')
    
    figure(7)
    plot(SS_VFId.gridk,cumsum(sum(SS_VFId.dist)),SS_VFIc.gridk,cumsum(sum(SS_VFIc.dist,1)),SS_PFI.gridk,cumsum(sum(SS_PFI.dist,1)),SS_EGM.gridk,cumsum(sum(SS_EGM.dist,1)),SS_CHEB.gridk,cumsum(sum(SS_CHEB.dist,1)),'LineWidth',2)  
    xlabel('assets','FontSize',14)
    ylabel('CDF','FontSize',14)
    title('asset distribution','FontSize',14)
    legend({'VFId','VFIc','PFI','EGM','CHEB'},'FontSize',12,'Location','SouthEast')

    figure(8)
    plot(SS_VFId.gridk,SS_VFId.V(1,:),SS_VFIc.gridk,SS_VFIc.V(1,:),SS_CHEB.gridk,SS_CHEB.V(1,:),'LineWidth',2)  
    xlabel('assets','FontSize',14)
    ylabel('value function','FontSize',14)
    title('value functions z_1','FontSize',14)
    legend({'VFId','VFIc','CHEB'},'FontSize',12,'Location','SouthEast')
    