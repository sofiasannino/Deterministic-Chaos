%%% DESCRIPTION -----------------------------------------------------------
%   the minimal vector of state variables defining a 1D periodic field


%%% INPUTS ----------------------------------------------------------------
%   u       the field in physical state (column vector of real numbers)
%   N       spatial resolution


%%% OUTPUTS ---------------------------------------------------------------
%   v       state vector (column vector of real numbers)


%%% REMARKS ---------------------------------------------------------------
%   1-  In forward integration of the KSE, the mean remains constant.
%       We assume the mean value is always zero, hence not stored.
%   2-  We assume that we work in the anti-symmetric subspace.
%       Thus, the Fourier coefficients are purely imaginary.


function v = field2vector(u,N)
    Nd = round(N/3);
    
    U = fft(u);
    v = imag(U(2:Nd+1));

end