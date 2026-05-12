%%% DESCRIPTION -----------------------------------------------------------
%   shooting-based Newton iterations for computing periodic orbits (POs)


%%% INPUTS ----------------------------------------------------------------
%   v0      guessed periodic point (column state vector)
%   T       guessed period of the PO
%   dt      reference time step for the required time marchings
%   nu      hyperviscosity parameter
%   N       spatial resolution


%%% OUTPUTS ---------------------------------------------------------------
%   v_best  exact periodic point (column state vector)
%   T_best  exact period of the PO
%   flag    flag == 1: search successful, otherwise: search failed

function [v_best,T_best,flag] = search4PO(v0,T,dt,nu,N)
    [~,k] = domain(N);
    U0 = fft(vector2field(v0,N));
    
    dUdt0 = -fft(ifft(U0,'symmetric').*ifft(complex(0,k).*U0,'symmetric'));
    dUdt0 = dealiase(dUdt0 + (k.^2 - nu * k.^4).*U0);
    dUdt0 = complex(0,imag(dUdt0));
    
    dvdt0 = field2vector(ifft(dUdt0,'symmetric'),N);
    
    options = optimoptions('fsolve','Display','iter','MaxIterations',75);
    [S,~,flag,~] = fsolve(@(X) recurrence(X),[v0;T],options);
    
    v_best = S(1:end-1);
    T_best = S(end);
    
    function F = recurrence(X)
        v0_ = X(1:end-1);
        T_ = X(end);
        
        [vT_,~] = KSE_integrate(v0_,T_,dt,0,nu,N);
        F = [vT_-v0_; dot(v0_-v0,dvdt0)];
    end
end