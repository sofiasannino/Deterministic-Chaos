function [u, U] = load_data(filename, Nsamples)
    data = read_file(strcat("../data/",filename));
    N = min(Nsamples, length(data));  % Only keep Nsamples data points
    data = data(1:N);
    U = mean(data);
    u = data - U;
end
