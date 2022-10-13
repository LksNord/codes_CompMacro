function [distp] = stepDist(dist,kpol,param,grid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to solve one step the distribution of HHs with continuous policies
% for reference see Young (2010, JEDC)
% inputs:
%       - param: structure of parameters
%       - grid: structure containing grids
%       - kopt: indices of optimal capital policy
%       - dist: current distribution
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

% initialize future distribution
    distp = zeros(param.nz,param.nkap);

% move distribution one step forward
    for i = 1:param.nkap      
        for j = 1:param.nz
             % for any combination of productivity and capital today
             % what is the distribution over z/k tomorrow 
             distp(:,indbelow(j,i)) = distp(:,indbelow(j,i)) + dist(j,i)*wbelow(j,i)*grid.Pz(j,:)';  
             distp(:,indabove(j,i)) = distp(:,indabove(j,i)) + dist(j,i)*wabove(j,i)*grid.Pz(j,:)';  
        end
    end
   

   