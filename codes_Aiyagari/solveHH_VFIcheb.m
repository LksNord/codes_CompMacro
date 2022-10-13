function [cpol,kpol,V] = solveHH_VFIcheb(r,w,param,grid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to solve the HH problem by VFI with projection method
% inputs:
%       - param: structure containing the necessary parameter values
%       - grid: structure containing grids
%       - r: interest rate
%       - w: wage rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define collocation points on [-1,1]
    xcolloc = NaN(1,param.Mcheb);
    for i=1:param.Mcheb
        xcolloc(i) = cos( pi/2 * ((2*i)-1)/param.Mcheb );
    end

    pssi = @(k) (k-grid.k(end))/(grid.k(end)-grid.k(1)); % function to map asset levels to [-1,1]
    
    kcolloc = grid.k(end) + (grid.k(end)-grid.k(1))*xcolloc; % asset levels associated with collocation points

    for i=1:param.Mcheb
        H(:,i) = cheby_vec(xcolloc(i),param.Ncheb); % chebyshev polynomials at the collocation nodes
    end
    
% define utility function
    if param.gamma==1
        u = @(c) log(c); % limit of CRRA for gamma=1
    else
        u = @(c) (c.^(1-param.gamma)-1)/(1-param.gamma); % CRRA
    end
 
% initialization
    err  = 1;                                % error
        
    weights0 = zeros(param.nz, param.Ncheb+1); % initial guess for chebychev weights
    weights1 = zeros(param.nz, param.Ncheb+1); % preallocation for updated chebychev weights
    kpolcheb = zeros(param.nz, param.Mcheb); % preallocation for asset choice
    V0 = weights0*H;                         % initial guess for value function
    V1 = zeros(param.nz, param.Mcheb);       % preallocation for updated value function

    while err > param.tol_vf
        
        for j=1:param.nz                                   % loop over idiosyncratic productivity
            for i=1:param.Mcheb
                Vint = @(kp) -( u(grid.z(j)*w + (1+r)*kcolloc(i)-kp) + param.beta*grid.Pz(j,:)*weights0*cheby_vec(pssi(kp),param.Ncheb)); % 'continuous' value function
                [kpolcheb(j,i),V1(j,i)] = fminbnd(Vint,grid.k(1),grid.z(j)*w + (1+r)*kcolloc(i)); % pick the optimal asset choice and store the associated value and policy
            end
        end
        
        V1=-V1; % invert from minimization problem
        
        for j=1:param.nz 
            weights1(j,:) = (H*H')\H*V1(j,:)';         % update weights
        end
        err = max(max(abs(V1-V0)));         % compute error
        V0 = V1;                            % update guess for value function
        weights0 = weights1;                % update guess for weights

    end

% find policy and value functions over actual grids
    kpol = zeros(param.nz, param.nkap); 
    V = zeros(param.nz, param.nkap);     
    for j=1:param.nz                                   
        for i=1:param.nkap
            Vint = @(kp) -( u(grid.z(j)*w + (1+r)*grid.k(i)-kp) + param.beta*grid.Pz(j,:)*weights0*cheby_vec(pssi(kp),param.Ncheb)); % 'continuous' value function
            [kpol(j,i),V(j,i)] = fminbnd(Vint,grid.k(1),grid.z(j)*w + (1+r)*grid.k(i)); % pick the optimal asset choice and store the associated value and policy
        end
    end
    V=-V;
    cpol = grid.z'*w + (1+r)*grid.k - kpol;