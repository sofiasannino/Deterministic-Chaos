%%% DESCRIPTION -----------------------------------------------------------
%   1D periodic field corresponding to a minimal vector of state variables


%%% INPUTS ----------------------------------------------------------------
%   v       state vector (column vector of real numbers)
%   N       spatial resolution


%%% OUTPUTS ---------------------------------------------------------------
%   u       the field in physical state (column vector of real numbers)


%%% REMARKS ---------------------------------------------------------------
%   1-  In forward integration of the KSE, the mean remains constant.
%       We assume the mean value is always zero, hence the mean value of
%       'u' is zero by construction.
%   2-  We do not fill the padded elements in Fourier space, hence no
%       additional dealiasing is needed.
%   3-  We assume that we work in the anti-symmetric subspace.
%       Thus, the Fourier coefficients are purely imaginary.


function u = vector2field(v,N)
    Nd = round(N/3);

    U = zeros(N,1);
    U(2:Nd+1) = v(1:Nd)*1j;    
    U(N-Nd+1:N) = flip(conj(U(2:Nd+1)));
    
    u = ifft(U,'symmetric');
end