function [cpol,kpol] = solveHH_PFI(r,w,param,grid,cguess)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to solve the HH problem by policy function iteration (PFI)
% inputs:
%       - param: structure containing the necessary parameter values
%       - grid: structure containing grids
%       - r: interest rate
%       - w: wage rate
%       - cguess: either 'false' or a matrix with appropriate sizes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% initialization
    err  = 1;                                % error
    if (cguess==false)
        cpol = (grid.z'*w + r*grid.k);       % guess for consumption policy
    else
        cpol = cguess;
    end
    cpol1= zeros(param.nz,param.nkap);       % preallocation for update
    wealth = grid.z'*w + (1+r)*grid.k;       % grid for cash on hand
    klarge = ones(param.nz,1)*grid.k;        % repeated asset grid
    
    while err > param.tol_pol

        for j=1:param.nz         % current state
            for i=1:param.nkap   % current asset holdings

                ready=0;

                % check if agent is borrowing constrained
                if ((wealth(j,i)-grid.k(1)))^(-param.gamma) >= param.beta*(1+r)*grid.Pz(j,:)*((cpol(:,1)).^(-param.gamma)) % current MU too high even at borrowing constraint
                    ready       = 1;
                    cpol1(j,i)  = (wealth(j,i)-grid.k(1));
                end

                % check if agent is savings constrained
                if ((ready==0) && (wealth(j,i)>grid.k(end)) &&  (((wealth(j,i)-grid.k(end)))^(-param.gamma) <= param.beta*(1+r)*grid.Pz(j,:)*(cpol(:,end)).^(-param.gamma))) % current MU too low even at max savings
                    ready       = 1;
                    cpol1(j,i)     = (wealth(j,i)-grid.k(end));
                end

                % if still ready==0, then asset decision is interior and we continue using bisection
                if (ready==0)
                    c_min       = max(0,(wealth(j,i)-grid.k(end))); % lowest possible consumption: zero or whatever left if choosing max saving
                    c_max       = (wealth(j,i)-grid.k(1));          % highest possible consumption: consume all cash on hand and max out on borrowiing going forward

                    err=1;
                    while (abs(err)>param.tol_pol)
                        cc      = c_min+0.5*(c_max-c_min);

                        % error in the Euler equation for current guess on todays consumption cc
                        err = cc^(-param.gamma)- param.beta*(1+r)*grid.Pz(j,:)*(interpLN((wealth(j,i)-cc),grid.k,cpol)).^(-param.gamma);
 
                        % update bisection bounds
                        % (this might be much faster with more intelligent
                        % updates, e.g. golden search or Newton/Quasi-Newton method)
                        if      (abs(err) > param.tol_pol && err<0) % todays MU too low --> decrease consumption guess
                            c_max   = cc;
                        elseif  (abs(err) > param.tol_pol && err>0) % todays MU too high --> increase consumption guess
                            c_min   = cc;
                        end

                    end

                    % store consumption solution for this gridpoint
                    cpol1(j,i) = cc;

                end

            end
        end
        
        err = max(max(abs(cpol-cpol1)));         % compute error
        cpol = cpol1;                            % update guess for consumption policy
    end

% compute asset policy
    kpol = wealth - cpol; 
