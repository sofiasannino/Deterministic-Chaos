function G = fit_gaussian(du, samples)
    % du      : vector of velocity increments [m]
    % samples : samples sampled according to prob. density function (PDF)
    % G       : A gaussian distribution fit G(du)

    mu = mean(samples);
    sigma = std(samples,1);

    G = 1 / (sqrt(2*pi)*sigma) * exp(-0.5 * ((du-mu) / sigma).^2);
end