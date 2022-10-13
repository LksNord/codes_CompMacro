function [dist,trans] = getDist_continuous(param,grid,kpol,trans_only)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to solve for the distribution of HHs with continuous policies
% for reference see Young (2010, JEDC)
% inputs:
%       - grid: structure containing grids
%       - kopt: indices of optimal capital policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% find index below capital policy
    indbelow = ones(param.nz,param.nkap);
    for i=2:param.nkap
        indbelow(kpol>=grid.k(i))=i;
    end
    indbelow = min(indbelow,param.nkap-1);
    indabove = indbelow + 1;
    
% find weights attached to point below/above policy
    wabove = (kpol-grid.k(indbelow))./(grid.k(indbelow+1)-grid.k(indbelow));
    wabove = min(wabove,1); wabove = max(wabove,0); % just to be safe, should not be binding
    wbelow= 1-wabove;

% fill transition matrix for the distribution
    %trans = zeros(param.nz*param.nkap,param.nz*param.nkap);  % preallocate generalized transition matrix
    trans =spalloc(param.nz*param.nkap,param.nz*param.nkap,param.nz*param.nkap*2*param.nz);
    % operating on a sparse matrix (one with many zeros) is MUCH faster in matlab!
    % try commenting in the alternative initialization of trans...
    
    for i = 1:param.nkap      
        for j = 1:param.nz
             % for any combination of productivity and capital today
             % (columns) what is the distribution over z/k tomorrow (rows)
             trans((i-1)*param.nz+j,(indbelow(j,i)-1)*param.nz+1:indbelow(j,i)*param.nz) = wbelow(j,i)*grid.Pz(j,:);  
             trans((i-1)*param.nz+j,(indabove(j,i)-1)*param.nz+1:(indabove(j,i))*param.nz) = wabove(j,i)*grid.Pz(j,:);  
        end
    end
    
    dist = NaN(param.nz,param.nkap);
    if (trans_only~=1)

        % Compute stationary distribution over productivities and assets
            probst = (1/(param.nz*param.nkap))*ones(1,param.nz*param.nkap);   % initialize stationary distribution (vector form)
            err = 1;                               % error for solving for the stationary distribution by iteration
            while err > 1e-14                
               probst1 = probst*trans;              % interate on stationary distribution with transition matrix
               err = max(abs(probst1-probst));      % compute error
               probst = probst1;                    % update guess
            end

        % reshape distribution to two dimensions
            dist=reshape(probst,param.nz,param.nkap);
    end  