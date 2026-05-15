%%% DESCRIPTION -----------------------------------------------------------
%   1D spatial domain in physical and Fourier states


%%% INPUTS ----------------------------------------------------------------
%   N   spatial resolution


%%% OUTPUTS ---------------------------------------------------------------
%   x   column vector of grid points
%   k   column vector of wave numbers


function [x,k] = domain(N)
    dx = 1./N;
    x = linspace(0,1.-dx,N)' - 1./2; %[0.5, 0.5)
    k = (2 * pi)*[0:N/2-1 0 -N/2+1:-1]';
end