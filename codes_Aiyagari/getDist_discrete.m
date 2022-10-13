function [dist] = getDist_discrete(param,grid,kopt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to solve for the distribution of HHs with discrete policies
% inputs:
%       - grid: structure containing grids
%       - kopt: indices of optimal capital policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% fill transition matrix for the distribution
    %trans = zeros(param.nz*param.nkap,param.nz*param.nkap);  % preallocate generalized transition matrix
    trans =spalloc(param.nz*param.nkap,param.nz*param.nkap,param.nz*param.nkap*2*param.nz);
    % operating on a sparse matrix (one with many zeros) is MUCH faster in matlab!
    % try commenting in the alternative initialization of trans...
    
    for i = 1:param.nkap      
        for j = 1:param.nz
             % for any combination of productivity and capital today
             % (columns) what is the distribution over z/k tomorrow (rows)
             trans((i-1)*param.nz+j,(kopt(j,i)-1)*param.nz+1:kopt(j,i)*param.nz) = grid.Pz(j,:);  
        end
    end
   
% Compute stationary distribution over productivities and assets
    probst = (1/(param.nz*param.nkap))*ones(1,param.nz*param.nkap);   % initialize stationary distribution (vector form)
    err = 1;                               % error for solving for the stationary distribution by iteration
    while err > 1e-10                
       probst1 = probst*trans;              % interate on stationary distribution with transition matrix
       err = max(abs(probst1-probst));      % compute error
       probst = probst1;                    % update guess
    end

% reshape distribution to two dimensions
    dist=reshape(probst,param.nz,param.nkap);
   