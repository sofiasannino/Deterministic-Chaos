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
nu = 6.5e-4;                    % hyperviscosity parameter
N = 128;                        % spatial resolution
t_study = 2500.0 * nu;          % analysis time period
dt = 0.1 * nu;                  % time step size for time integration
dt_store = 1.0 * nu;            % time intervals of storing a snapshot
T_max = 120.0 * nu;             % maximum T for recurrent flow analysis
epsilon = 1e-6;                 % perturbation amplitude for computing J

%% Initial condition
[x,~] = domain(N);              % construct the spatial domain
u0 = sin(2.0*pi*x) / sqrt(nu);  % initial condition in physical state
v0 = field2vector(u0,N);        % initial state vector

%% Time integrate the KSE
[V,t_vec] = KSE_integrate(v0,t_study,dt,dt_store,nu,N);

%% Recurrence indicator
% Maximum lag in number of stored snapshots
max_lag = floor(T_max / dt_store);

% Make sure we do not exceed the available trajectory length
max_lag = min(max_lag, length(t_vec)-1);

% Time delays T = j * dt_store, with T > 0
T_vec = (1:max_lag)' * dt_store;

% We can only use times t such that t + T_max is still inside the data
n_t = length(t_vec) - max_lag;
t_rec = t_vec(1:n_t);

% Allocate recurrence matrix
R = zeros(n_t, max_lag);

for j = 1:max_lag
    % Difference v(t+T) - v(t) for all admissible t
    dV = V(:, (1+j):(n_t+j)) - V(:, 1:n_t);

    % Norms column by column
    numerator = vecnorm(dV, 2, 1);
    denominator = vecnorm(V(:, 1:n_t), 2, 1);

    R(:, j) = numerator ./ denominator;
end

% Avoid log(0), only for numerical safety
R(R <= 0) = eps;

% Plot ln(r(t,T))
plot2G = figure('Name', 'Plot2G');
imagesc(t_rec/nu, T_vec/nu, log(R)');
set(gca, 'YDir', 'normal');
colorbar;
colormap(parula);

xlabel('$t/\nu\;\mathrm{[s^2/m^2]}$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

ylabel('$T/\nu\;\mathrm{[s^2/m^2]}$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

% Save the plot automatically
if isgraphics(plot2G)
    exportgraphics(plot2G, fullfile(params.plot.folder, "Plot2G.png"), "Resolution", 300);
else
    disp('Skipped exporting plot2G (figure was closed).');
end

%% Manually selected recurrence pairs
% Hardcode here the selected pairs from Plot G
% Columns: [t/nu, T/nu]
pairs_scaled = [   
    1552   77;
    813     52;
    1435   57;
    954    106;
    1354   25
];

% Convert back to physical time
t_guess = pairs_scaled(:,1) * nu;
T_guess = pairs_scaled(:,2) * nu;

% Snap to nearest stored snapshots
idx_guess = zeros(size(t_guess));
T_guess_grid = zeros(size(T_guess));

for q = 1:length(t_guess)
    [~, idx_guess(q)] = min(abs(t_vec - t_guess(q)));
    [~, iT] = min(abs(T_vec - T_guess(q)));

    T_guess_grid(q) = T_vec(iT);
end

% Corresponding initial guess states
v_guess = V(:, idx_guess);

n_guesses = size(pairs_scaled, 1);

for q = 1:size(pairs_scaled,1)
    [~, it] = min(abs(t_rec - pairs_scaled(q,1)*nu));
    [~, iT] = min(abs(T_vec - pairs_scaled(q,2)*nu));

    r_selected(q,1) = R(it,iT);
    log_r_selected(q,1) = log(R(it,iT));
end

%% Converge UPO guesses
v_upo = zeros(size(V,1), n_guesses);
T_upo = zeros(n_guesses, 1);
flag_upo = zeros(n_guesses, 1);
closure_error_upo = zeros(n_guesses, 1);

for q = 1:n_guesses
    % Initial guess from recurrence plot
    vg = v_guess(:, q);
    Tg = T_guess_grid(q);

    % Converge periodic orbit
    [v_best, T_best, flag] = search4PO(vg, Tg, dt, nu, N);

    % Store result
    v_upo(:, q) = v_best;
    T_upo(q) = T_best;
    flag_upo(q) = flag;

    % Check closure error
    v_end = KSE_integrate(v_best, T_best, dt, 0, nu, N);
    closure_error_upo(q) = norm(v_end - v_best) / norm(v_best);
end


%% Plot H: chaotic attractor with UPOs
% Project chaotic trajectory
[E_att, P_att, D_att] = projection(V, N, nu);

plot2H = figure('Name', 'Plot2H');
hold on; grid on;

% Plot chaotic attractor
plot3(E_att, P_att, D_att, ...
    'Color', [0.75 0.75 0.75], ...
    'LineWidth', 0.8, ...
    'DisplayName', 'Chaotic attractor');

% Plot each computed UPO
for q = 1:n_guesses
    % Skip candidates with invalid period/state
    if any(isnan(v_upo(:,q))) || isnan(T_upo(q)) || T_upo(q) <= 0
        continue
    end

    % Time-march the converged UPO for one full period
    [V_orbit, ~] = KSE_integrate(v_upo(:,q), T_upo(q), dt, dt_store, nu, N);

    % Project UPO onto (E,P,D)
    [E_orb, P_orb, D_orb] = projection(V_orbit, N, nu);

    % Overlay UPO
    plot3(E_orb, P_orb, D_orb, ...
        'LineWidth', 1.8, ...
        'DisplayName', ['$\mathrm{UPO}\ ', num2str(q), '$']);
end

xlabel('$E$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

ylabel('$P$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

zlabel('$D$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

legend('Interpreter', 'latex', 'Location', 'best', ...
    'FontSize', params.plot.legendFontSize);

view(3);

% Save the plot automatically
if isgraphics(plot2H)
    exportgraphics(plot2H, fullfile(params.plot.folder, "Plot2H.png"), "Resolution", 300);
else
    disp('Skipped exporting plot2H (figure was closed).');
end

%% Floquet multipliers and time scales 
lambda_max_upo = zeros(n_guesses, 1);
tF_star = zeros(n_guesses, 1);
stability_upo = strings(n_guesses, 1);

for q = 1:n_guesses
    % Compute Jacobian of the flow map over one UPO period
    J = Jacobian(v_upo(:, q), T_upo(q), epsilon, dt, nu, N);

    % Floquet multipliers
    lambda = eig(J);
    lambda_abs = abs(lambda);

    % Ignore the neutral multiplier closest to 1
    [~, idx_neutral] = min(abs(lambda - 1));
    lambda_abs(idx_neutral) = [];

    % Leading non-neutral Floquet multiplier
    lambda_max_upo(q) = max(lambda_abs);

    % Stability classification
    if lambda_max_upo(q) > 1
        stability_upo(q) = "unstable";
    else
        stability_upo(q) = "stable";
    end

    % Floquet divergence time scale
    tF_star(q) = T_upo(q) / log(lambda_max_upo(q));
end

%% --- Display results ---
fprintf('====================================================== Part 2.4 results ======================================================\n');
fprintf('%-5s  %-12s  %-12s  %-13s  %-6s  %-12s  %-12s  %-12s  %-10s  %-12s\n', ...
    'UPO', 't', 'T_guess', 'log r(t,T)', 'Flag', 'T_UPO', 'Closure Err', 'Max |lambda|', 'Stability', 'tF_star');

for q = 1:n_guesses
    fprintf('%-5d  %-12e  %-12.4e  %-12e  %-6d  %-12e  %-12e  %-12e  %-10s  %-12e\n', ...
        q, pairs_scaled(q,1)*nu, T_guess_grid(q), log_r_selected(q), flag_upo(q), ...
        T_upo(q), closure_error_upo(q), lambda_max_upo(q), stability_upo(q), tF_star(q));
end