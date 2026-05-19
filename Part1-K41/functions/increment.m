function [x, d] = increment(u, l, f, U)
    % u: vector of perturbative velocity [m/s]
    % l: length for the increment calculation [m]
    % f: frequency of data sampling [Hz]
    % U: mean flow speed [m/s]
    % x: vector of lengths [m]
    % d: vector of velocity increments [m/s]
    
    N = length(u);

    % Consistently with Part 1, use Taylor's hypothesis as x = - U t.
    % Therefore, increasing the time index corresponds to more
    % negative x values
    dx = abs(U) / f;
    lag = round(l/dx);

    % With x_i = -(i-1)dx, the point x+l corresponds to an earlier index
    d = u(1:end-lag) - u(1+lag:end);
    x = -((lag):(N-1)) * dx;
end