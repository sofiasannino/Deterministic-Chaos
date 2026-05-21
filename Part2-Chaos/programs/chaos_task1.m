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
nu_values = [1e-2, 1e-3, 8 * 1e-4, 6.5 * 1e-4, 2 * 1e-4];
nu = nu_values(1);             % hyperviscosity parameter

N = 128;                       % spatial resolution
T_trans = 1000 * nu;           % transient time period
T_int = 250 * nu;              % analysis time period
dt = 0.1 * nu;                 % time step size for time integration
dt_store = 0.1 * nu;           % time intervals of storing a snapshot
epsilon = 1e-2;                % relative amplitude of perturbation

%% Step 1: initial condition
[x,~] = domain(N);             % construct the spatial domain
u0 = sin(2.0*pi*x) / sqrt(nu); % initial condition in physical state
v0 = field2vector(u0,N);       % initial state vector

%% Step 2: transient time integration
[v1000,~] = KSE_integrate(v0,T_trans,dt,0,nu,N);

%% Step 3: perturbing the state vector
r = epsilon * (2 * rand(size(v1000)) - 1) .* v1000;
v_perturbed = v1000 + r;       % Add perturbation to the starting state vector

%% Steps 2 & 3: Advance both states for T_int
[vv1, tt1] = KSE_integrate(v1000, T_int, dt, dt_store, nu, N);
[vv2, tt2] = KSE_integrate(v_perturbed, T_int, dt, dt_store, nu, N);

% Convert state vector to physical fields for plotting
u1 = zeros(N, length(tt1));
u2 = zeros(N, length(tt2));

for i = 1:length(tt1)
    u1(:,i) = vector2field(vv1(:,i), N);
    u2(:,i) = vector2field(vv2(:,i), N);
end

%% Steps 2, 3, 4, 5: Plotting A(nu)
plotname = "Plot2A";
extraplotname = "Plot2Extra";

for i = 1:length(nu_values)
    if nu == nu_values(i)
        plotname = plotname + num2str(i) + ".png";
        extraplotname = extraplotname + num2str(i) + ".png";
        break;
    end
end

plot2A = figure("Name", plotname);

% Panel 1: Space-time contour of u1(x,t)
subplot(2, 2, 1);
contourf(x, tt1/nu, u1', 50, 'LineStyle', 'none');
colorbar;

xlabel('$x\;[\mathrm{m}]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

ylabel('$\mathrm{t / \nu}\;[\mathrm{s^2/m^2}]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

title('Unperturbed u_1(x,t)', ...
    'FontSize', params.plot.fontSize);

% Panel 2: Space-time contour of u2(x,t)
subplot(2, 2, 2);
contourf(x, tt2/nu, u2', 50, 'LineStyle', 'none');
colorbar;

xlabel('$x\;[\mathrm{m}]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

ylabel('$\mathrm{t / \nu}\;[\mathrm{s^2/m^2}]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

title('Perturbed u_2(x,t)', ...
    'FontSize', params.plot.fontSize);

% Panel 3: Absolute difference |u1(x,t) - u2(x,t)|
subplot(2, 2, 3);
abs_diff = abs(u1 - u2);
contourf(x, tt1/nu, abs_diff', 50, 'LineStyle', 'none');
colorbar;

xlabel('$x\;[\mathrm{m}]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

ylabel('$\mathrm{t / \nu}\;[\mathrm{s^2/m^2}]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

title('Difference |u_1 - u_2|', ...
    'FontSize', params.plot.fontSize);

% Panel 4: 2D Projection onto P and D
[E1, P1, D1] = projection(vv1, N, nu);
[E2, P2, D2] = projection(vv2, N, nu);

subplot(2, 2, 4);
hold on; grid on;

color_red_alpha  = [1, 0, 0, 0.4]; 
color_blue_alpha = [0, 0, 1, 0.4];

plot(P2, D2, 'Color', color_red_alpha, 'LineWidth', 1.5, 'DisplayName', 'Perturbed');
plot(P1, D1, 'Color', color_blue_alpha, 'LineWidth', 1.5, 'DisplayName', 'Unperturbed');
plot(P1(1), D1(1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'DisplayName', 'Start Point');

xlabel('Energy Production Rate (P)', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

ylabel('Energy Dissipation Rate (D)', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

title('2D Phase Portrait', ...
    'FontSize', params.plot.fontSize);

legend('Location', 'best');

%% Extra: compute divergence timescale
error_norm = vecnorm(vv1 - vv2);

% Plot the error
plot2Extra1 = figure('Name', extraplotname);
ax_extra = gca;

semilogy(ax_extra, tt1/nu, error_norm, 'k', 'LineWidth', 1.5);
hold(ax_extra, 'on'); grid(ax_extra, 'on');

xlabel(ax_extra, '$\mathrm{t / \nu}\;[\mathrm{s^2/m^2}]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

ylabel(ax_extra, '$|| v_1(t) - v_2(t) ||\;[\mathrm{m/s}]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

title(ax_extra, 'Divergence of perturbed trajectory', ...
    'FontSize', params.plot.fontSize);

mean_state_size = mean(vecnorm(vv1));
threshold = 0.1 * mean_state_size;
threshold_label = '10% Mean Attractor Size';

% Alternative
% threshold = 100 * error_norm(1);
% threshold_label = '100x Initial Error';

yline(ax_extra, threshold, 'r--', threshold_label, 'LabelHorizontalAlignment', 'left');
dev_idx = find(error_norm >= threshold, 1);
if ~isempty(dev_idx)
    t_deviation = tt1(dev_idx) / nu;
    plot(ax_extra, t_deviation, error_norm(dev_idx), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end

% Save the plots automatically
if isgraphics(plot2A)
    exportgraphics(plot2A, fullfile(params.plot.folder, plotname), "Resolution", 300);
else
    disp('Skipped exporting ' + plotname + ' (figure was closed).');
end

if isgraphics(plot2Extra1)
    exportgraphics(plot2Extra1, fullfile(params.plot.folder, extraplotname), "Resolution", 300);
else
    disp('Skipped exporting ' + extraplotname + ' (figure was closed).');
end

%% --- Display results ---
fprintf('================ Part 2.1 results ================\n');
if ~isempty(dev_idx)
    t_deviation = tt1(dev_idx) / nu;
    fprintf('For nu = %g, significant deviation occurs at t/nu approx: %e\n', nu, t_deviation);
else
    fprintf('For nu = %g, no significant deviation occurred (system is stable).\n', nu);
end