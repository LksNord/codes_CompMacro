function [Finterp] = interpLN(x,xgrid,F)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolation routine
% inputs:
%       - x: point to interpolate to
%       - xgrid: grid for input function
%       - F: function values over grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This user written function is simpler and hence faster than the matlab
% buildin, it is sufficient for the simple problem at hand. If you are
% looking for something more robust (but slower) check e.g. Matlabs buildin
% function "interp1()".

    xlow = max(sum(x>xgrid),1);        % grid point just below choice for k
    xhigh = xlow + 1;                   % grid point just above choice for k
    
    % linear interpolation of current value function guess
    Finterp = F(:,xlow) + (x-xgrid(xlow))*(F(:,xhigh) - F(:,xlow))/(xgrid(xhigh) - xgrid(xlow));

    
