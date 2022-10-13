function [cpol,kpol,V0] = solveHH_VFIc(r,w,param,grid,Vguess)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to solve the HH problem by discrete VFI
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
 
% initialization
    err  = 1;                                % error
    if (Vguess==false)  
        V0 = u(grid.z'*w + r*grid.k)/(1-param.beta);  % guess for value function
    else
        V0 = Vguess;
    end
    
    V1   = zeros(param.nz,param.nkap);       % preallocation for updated value function
    kpol = zeros(param.nz,param.nkap);       % preallocation for asset choice
    wealth = grid.z'*w + (1+r)*grid.k;       % grid for cash on hand

    while err > param.tol_vf
        
        for j=1:param.nz                                   % loop over idiosyncratic productivity
            for i=1:param.nkap
                Vint = @(kp) -( u(wealth(j,i)-kp) + param.beta*grid.Pz(j,:)*interpLN(kp,grid.k,V0)); % 'continuous' value function
                [kpol(j,i),V1(j,i)] = fminbnd(Vint,grid.k(1),min(wealth(j,i),grid.k(end))); % pick the optimal asset choice and store the associated value and policy
            end
        end
        
        V1=-V1; % invert from minimization problem
        
        err = max(max(abs(V1-V0)));         % compute error
        V0 = V1;                            % update guess for value function
    end

cpol = grid.z'*w + (1+r)*grid.k - kpol;