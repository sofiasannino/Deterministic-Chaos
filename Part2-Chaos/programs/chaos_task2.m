%% --- Move to this file's folder ---
try
    cd(fileparts(matlab.desktop.editor.getActiveFilename));
catch
    cd(fileparts(mfilename('fullpath')));
end

%% --- Clear workspace ---
clc;
clear;
addpath("../functions/");

%% --- Output folder for plots ---
params.plot.folder = "../plots";

if ~exist(params.plot.folder, "dir")
    mkdir(params.plot.folder);
end

%% inputs
nu_values = [8 * 1e-4, 6.5 * 1e-4, 2 * 1e-4]; %viscosity values 
nu = nu_values(3);             % hyperviscosity parameter
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
log_sum = zeros(N_exp, 1);     % vector to keep the logarithmic sums 

% Initialize the plot variables
plotA = figure();
axA = axes(plotA);

grid(axA, 'on'); 
hold(axA, 'on');
xlabel(axA, 't/\nu'); 
ylabel(axA, '$\tilde{\chi_i}$', 'Interpreter', 'latex');
line_colors = lines(N_exp); 
ani_lines = gobjects(N_exp, 1);

for q = 1:N_exp
    labelStr = ['$\tilde{\chi}_' num2str(q) '$'];
    ani_lines(q) = animatedline(axA, 'LineWidth', 2, 'Color', line_colors(q,:), ...
                                'DisplayName', labelStr);
end

legend(axA, 'Interpreter', 'latex', 'Location', 'best');

for i = 1:N_norm
    J = Jacobian(v0,tau,epsilon,dt,nu,N);
    
    % Map the deviation vectors
    V = J * Q;

    % Compute QR decomposition V = QR, update orthonormal matrix
    [Q, R] = qr(V , "econ");
    
    % Update
    log_sum = log_sum + log(abs(diag(R)));

    [v0,~] = KSE_integrate(v0,tau,dt,0,nu,N);
    t(i) = i*tau;
    X(:, i) = log_sum./t(i);

    for q = 1:N_exp
        addpoints(ani_lines(q), t(i)/nu, nu.*X(q,i));
    end
    
    % Render the updates to the screen every 50 steps
    if rem(i, 50) == 0
        drawnow limitrate;
    end
end

% Save the plot automatically
if isgraphics(plotA)
    name = "Chaos2A" + nu + ".png";
    exportgraphics(plotA, fullfile(params.plot.folder, name), "Resolution", 300);
else
    disp('Skipped exporting plotA (figure was closed).');
end

% Print results
fprintf('\nFinal rescaled Lyapunov exponents (nu * chi_i) at t/nu = %.2f:\n', t(end)/nu);
final_rescaled_X = nu .* X(:, end);
for q = 1:N_exp
    fprintf('  chi_%d = %g\n', q, final_rescaled_X(q));
end

fprintf('\nt^*_L = %.2f:\n', 1/max(X(:, end)));