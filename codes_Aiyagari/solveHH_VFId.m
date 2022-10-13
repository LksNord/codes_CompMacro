function [cpol,kpol,kopt,V0] = solveHH_VFId(r,w,param,grid,Vguess)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to solve the HH problem by continuous VFI
% inputs:
%       - param: structure containing the necessary parameter values
%       - grid: structure containing grids
%       - r: interest rate
%       - w: wage rate
%       - Vguess: either 'false' or a matrix with appropriate sizes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define utility function
if param.gamma==1
    u = @(c) log(c); % limit of CRRA for gamma=1
else
    u = @(c) (c.^(1-param.gamma)-1)/(1-param.gamma); % CRRA
end

% Fill matrix of flow utilities
    utilm = zeros(param.nz,param.nkap,param.nkap);              % preallocate matrix to store instantaneous utilities
    
    for j=1:param.nz          % iterate over idiosyncratic productivity
            cons = grid.z(j)*w + (1+r)*grid.k' -grid.k;     % cons is vector of consumption levels given current state and all possible asset choices
            util = u(cons);                                     % apply utility function
            util(cons <= 0) = -1e8;                             % rule out negative consumption - allocate large negative number --> never chosen        
            utilm(j,:,:)=util;                                  % store in matrix
    end

% initialization
    err  = 1;                                % error
    if (Vguess==false)
        V0 = u(grid.z'*w + r*grid.k)/(1-param.beta);  % guess for value function   
    else
        V0 = Vguess;
    end
    
    V1   = zeros(param.nz,param.nkap);       % preallocation for updated value function
    kopt = zeros(param.nz,param.nkap);       % preallocation for index of asset choice

    while err > param.tol_vf
        
        for j=1:param.nz                                   % loop over idiosyncratic productivity

                util = squeeze(utilm(j,:,:));               % instantaneous utilities for current state and all possible choices            
                EV = grid.Pz(j,:)*V0;                       % expected future value function
                Vint = util + param.beta*EV;                % value function for all possible choices          
                [V1(j,:),kopt(j,:)] = max(Vint,[],2);          % pick the optimal asset choice and store the associated value and policy

        end
        
        err = max(max(abs(V1-V0)));         % compute error
        V0 = V1;                            % update guess for value function
    end

kpol = grid.k(kopt);  % transfor policy indices into values
cpol = grid.z'*w + (1+r)*grid.k - kpol;