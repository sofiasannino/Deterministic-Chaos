function S = structure_function(u, lvec, f, U, n)
    % u    : vector of perturbative velocity [m/s]
    % lvec : length for the increment calculation [m]
    % f    : frequency of data sampling [Hz]
    % U    : mean flow speed [m/s]
    % n    : order of the structure function [-]
    % S    : vector of the structure function [m^n/s^n]

    % We have a time signal taken at a fixed point r_probe. It is
    % interpreted as a spatial signal through Taylor's
    % frozen-flow hypothesis using x = -U t.
    % Therefore, increasing the time index corresponds
    % to moving toward smaller x.
    
    % We furthermore assume isotropy of the structure function, meaning that
    % the difference in velocities depends only on the modulus of l. This
    % allows us to take the vectorial increment l in the direction of flow
    % (where we can use Taylor's hp. to estimate u(r_probe+l)).
    
    % The final assumption we need to make is that of ergodicity. We
    % calculate what should rigorously be statistical means as time
    % averages instead.
    
    % With this convention, u(j+Delta j) corresponds to u(x-l).
    % delta u_parallel(x,l) = u(x+l) - u(x),
    % therefore
    % delta u_parallel = u(j) - u(j+Delta j).

    u = u(:);
    lvec = lvec(:);

    S = NaN(length(lvec), 1);

    % spatial distance between consecutive samples
    dx = abs(U) / f;

    for i = 1:length(lvec)

        index_shift = round(lvec(i) / dx);

        if index_shift < 1
            warning('l = %.3e m is smaller than the maximum resolution possible dx = %.3e m.', lvec(i), dx);
            continue;
        end

        if index_shift >= length(u)
            warning('l = %.3e m is too large for the available signal length.', lvec(i));
            continue;
        end

        % u(1:end-index_shift) corresponds to u(x+l)
        % u(1+index_shift:end) corresponds to u(x)
        du_parallel = u(1:end-index_shift) - u(1+index_shift:end);

        S(i) = mean(du_parallel.^n);

    end
end

