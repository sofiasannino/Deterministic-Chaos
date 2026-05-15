close all; clear; clc

%% inputs
addpath('../functions/')       % folder containing functions
nu_values = [1e-2, 1e-3, 8 * 1e-4, 6.5 * 1e-4, 2 * 1e-4]; %viscosity values 
nu = nu_values(2);             % hyperviscosity parameter
N = 128;                       % spatial resolution
dt = 0.1 * nu;                 % time step size for time integration
epsilon = 1e-6;                % perturbation amplitude for computing J
tau = 2 * nu;                  % time intervals of computing J
N_norm = 5000;                 % number of re-normalizations
N_exp = 5;                     % number of Lyapunov exponents

%% initial condition
[x,~] = domain(N);             % construct the spatial domain
u0 = sin(2.0*pi*x) / sqrt(nu); % initial condition in physical state
v0 = field2vector(u0,N);       % initial state vector
n = length(v0);

%% compute Lyapunov exponents
Q = zeros(length(v0),N_exp);   % allocate matrix of orthonormal vectors
for i = 1:N_exp
    Q(i,i) = 1;                % initialize orthonormal vectors
end

X = zeros(N_exp,N_norm);       % allocate matrix of Lyapunov exp. history
t = zeros(N_norm,1);           % allocate vector of normalization instances
log_sum = zeros(N_exp, 1); % vector to keep the logarithmic sums 
figure
for i = 1:N_norm
    J = Jacobian(v0,tau,epsilon,dt,nu,N);
    
    % map the deviation vectors
    V = J * Q; % n*N_exp

    % compute QR decomposition V = QR, update orthonormal matrix
    [Q, R] = qr(V , "econ");
    
    % update
    log_sum = log_sum + log(abs(diag(R)));

    [v0,~] = KSE_integrate(v0,tau,dt,0,nu,N);
    t(i) = i*tau;
    
    X(:, i) = log_sum./t(i);

    if rem(i,20) == 0          % update figure every 20 re-normalizations
        clf; grid on; hold on
        for q = 1:N_exp
            plot(t(1:i)/nu,nu.*X(q,1:i),'LineWidth',2)
        end
        xlabel('t/\nu'); ylabel('$\tilde{\chi_i}$', 'Interpreter', 'latex')
        drawnow
    end
end
