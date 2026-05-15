function [x, d] = increment(u, l, f, U)
    % u: vector of perturbative velocity [m/s]
    % l: length for the increment calculation [m]
    % f: frequency of data sampling [Hz]
    % U: mean flow speed [m/s]
    % x: vector of lengths [m]
    % d: vector of velocity increments [m/s]
    
    N = length(u);

    % Convention: consistently with Part 1, we use Taylor's hypothesis as
    % x = - U t. Therefore, increasing the time index corresponds to more
    % negative x values.
    dx = abs(U) / f;

    lag = round(l/dx);

    % With x_i = -(i-1)dx, the point x+l corresponds to an earlier index.
    % Therefore this implements delta u_parallel(x,l) = u(x+l) - u(x).
    d = u(1:end-lag) - u(1+lag:end);
<<<<<<< HEAD

    x = -((lag):(N-1)) * dx;

=======

    x = -((lag):(N-1)) * dx;
>>>>>>> 2891bd1ef2767442f22b0032f1d0c66c9aa9416e

end