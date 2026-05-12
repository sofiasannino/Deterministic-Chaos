%%% DESCRIPTION -----------------------------------------------------------
%   dealiasing a 1D periodic field in spectral state using the 2/3 rule


%%% INPUTS ----------------------------------------------------------------
%   U   vector of Fourier coefficients (vector of complex numbers)


%%% OUTPUTS ---------------------------------------------------------------
%   U_  the same as 'U' with 1/3 of the highest frequency modes set to 0


function U_ = dealiase(U)
    N = length(U);
    Nd = round(N/3);
    
    U_ = U;
    U_(Nd+2:N-Nd) = 0;
end