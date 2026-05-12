function [x, d] = increment(u, l, f, U)
    % u: vector of perturbative velocity [m/s]
    % l: length for the increment calculation [m]
    % f: frequency of data sampling [Hz]
    % U: mean flow speed [m/s]
    % x: vector of lengths [m]
    % d: vector of velocity increments [m/s]
    
    N = length(u);

    dx = U / f;

    lag = round(l/dx);

    d = u(1+lag:end) - u(1:end-lag);

    x = (0:length(d)-1)* dx;

end