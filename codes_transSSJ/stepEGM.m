function [cpol,kpol] = stepEGM(cpolp,r,rp,w,tau,T,param,grid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to solve one step of the HH problem by endogenous grid method (EGM)
% for reference see Carroll (2006, Economics Letters)
% inputs:
%       - param: structure containing the necessary parameter values
%       - grid: structure containing grids
%       - r: interest rate
%       - rp: next periods interest rate
%       - w: wage rate
%       - tau: labor income tax
%       - T: transfer
%       - cpolp: next periods consumption policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize asset policy
    kpol = zeros(param.nz,param.nkap);

% compute expected future marginal utility of consumption
    Ec = grid.Pz*cpolp.^(-param.gamma);

% backout current consumption
    c_impl = ((1+rp)*param.beta*Ec).^(-1/param.gamma);

% get current assets
    k_impl = ( c_impl + grid.k - (1-tau)*grid.z'*w - T)/(1+r);

% interpolate asset policy
    for j = 1:param.nz
        kpol(j,:) = interp1(k_impl(j,:),grid.k,grid.k,'linear','extrap'); 
        % I am using the buildin matlab function for clarity, you
        % might be able to code this yourself faster
    end

% correct for binding borrowing limit
    kpol(kpol<grid.k(1)) = grid.k(1);

% back out consumption from budget constraint
    cpol = (1-tau)*grid.z'*w + T + (1+r)*grid.k - kpol;


