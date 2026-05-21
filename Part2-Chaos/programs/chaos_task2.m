%% --- Clear workspace ---
clc;
clear;

%% --- Robust Path Setup ---
% Find the directory of this script, whether running via F5 or F9
if isempty(mfilename)
    % Fallback for when running code line-by-line (F9) in the Editor
    scriptDir = fileparts(matlab.desktop.editor.getActiveFilename);
else
    % Standard execution (F5 or command line)
    scriptDir = fileparts(mfilename('fullpath'));
end

% Build absolute paths relative to the script's directory
functionsDir = fullfile(scriptDir, '..', 'functions');
params.plot.folder = fullfile(scriptDir, '..', 'plots');

% Add dependencies to path
addpath(functionsDir);

% Create the plot folder if it doesn't exist
if ~exist(params.plot.folder, 'dir')
    mkdir(params.plot.folder);
end

%% --- Common setup ---
params.plot.fontSize = 20;
params.plot.legendFontSize = 18;
params.plot.insetFontSize = 16;
params.plot.lw = 2.1;
params.plot.lwsmall = 1.8;

%% Inputs
nu_values = [8 * 1e-4, 6.5 * 1e-4, 2 * 1e-4];
nu = nu_values(1);             % hyperviscosity parameter

N = 128;                       % spatial resolution
dt = 0.1 * nu;                 % time step size for time integration
epsilon = 1e-6;                % perturbation amplitude for computing J
tau = 2 * nu;                  % time intervals of computing J
N_norm = 5000;                 % number of re-normalizations
N_exp = 5;                     % number of Lyapunov exponents

%% Initial condition
[x,~] = domain(N);             % construct the spatial domain
u0 = sin(2.0*pi*x) / sqrt(nu); % initial condition in physical state
v0 = field2vector(u0,N);       % initial state vector
n = length(v0);

%% Compute Lyapunov exponents
Q = zeros(length(v0), N_exp);  % allocate matrix of orthonormal vectors
for i = 1:N_exp
    Q(i,i) = 1;                % initialize orthonormal vectors
end

X = zeros(N_exp, N_norm);      % allocate matrix of Lyapunov exp. history
t = zeros(N_norm, 1);          % allocate vector of normalization instances
log_sum = zeros(N_exp, 1);     % vector to keep the logarithmic sums 

% Initialize the plot variables
plotname = "Plot2B";

for i = 1:length(nu_values)
    if nu == nu_values(i)
        plotname = plotname + num2str(i) + ".png";
        break;
    end
end

plot2B = figure("Name", plotname);

axB = axes(plot2B);
grid(axB, 'on'); 
hold(axB, 'on');

xlabel('$\mathrm{t / \nu}\;[\mathrm{s^2/m^2}]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

ylabel('$\tilde{\chi_i}$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

line_colors = lines(N_exp); 
ani_lines = gobjects(N_exp, 1);

for q = 1:N_exp
    labelStr = ['$\tilde{\chi}_' num2str(q) '$'];
    ani_lines(q) = animatedline(axB, 'LineWidth', params.plot.lw, ...
                                'Color', line_colors(q,:), ...
                                'DisplayName', labelStr);
end

legend(axB, 'Interpreter', 'latex', ...
    'FontSize', params.plot.legendFontSize, ...
    'Location', 'northeast');

for i = 1:N_norm
    J = Jacobian(v0,tau,epsilon,dt,nu,N);
    
    % Map the deviation vectors
    V = J * Q;

    % Compute QR decomposition V = QR, update orthonormal matrix
    [Q, R] = qr(V , "econ");
    
    % Update the logarithmic sum
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
if isgraphics(plot2B)
    exportgraphics(plot2B, fullfile(params.plot.folder, plotname), "Resolution", 300);
else
    disp('Skipped exporting ' + plotname + ' (figure was closed).');
end

%% --- Display results ---
fprintf('================ Part 2.2 results ================\n');
fprintf('Final rescaled Lyapunov exponents (nu * chi_i) at t/nu = %e:\n', t(end)/nu);

final_rescaled_X = nu .* X(:, end);
for q = 1:N_exp
    fprintf('chi_%d = %e\n', q, final_rescaled_X(q));
end

fprintf('\nt^*_L = %e:\n', 1/max(X(:, end)));