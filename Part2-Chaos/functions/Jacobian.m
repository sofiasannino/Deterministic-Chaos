%%% DESCRIPTION -----------------------------------------------------------
%   Jacobian of the KSE flow map: J=dv(t)/dv(0)


%%% INPUTS ----------------------------------------------------------------
%   v0          reference point of the Jacobian (column state vector)
%   t           integration time interval
%   epsilon     perturbation magnitude for finite difference derivatrives
%   dt          step size in time integrations
%   nu          hyperviscosity parameter
%   N           spatial resolution


%%% OUTPUTS ---------------------------------------------------------------
%   J           Jacobian matrix

function J = Jacobian(v0,t,epsilon,dt,nu,N)

    n = length(v0);
    J = zeros(n,n);

    % central finite difference approximation of the Jacobian
    for j = 1:n

        % perturbation of the j-th component
        dv = zeros(n,1);
        dv(j) = epsilon;

        % perturbed initial conditions
        v_plus  = v0 + dv;
        v_minus = v0 - dv;

        % integrate the two perturbed trajectories
        [vT_plus,~]  = KSE_integrate(v_plus,t,dt,0,nu,N);
        [vT_minus,~] = KSE_integrate(v_minus,t,dt,0,nu,N);

        % j-th column of the Jacobian
        J(:,j) = (vT_plus - vT_minus) / (2 * epsilon);

    end

end

