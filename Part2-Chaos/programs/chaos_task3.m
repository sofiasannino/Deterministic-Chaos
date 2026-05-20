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

%% --- Common setup ---
params.plot.fontSize = 20;
params.plot.legendFontSize = 18;
params.plot.insetFontSize = 16;
params.plot.lw = 2.1;
params.plot.lwsmall = 1.8;

%% Inputs
N = 128;                       % spatial resolution
epsilon = 1e-6;                % perturbation amplitude for computing J
m = 20;                        % number of continuation steps
nu0 = 8e-4;                    % initial hyperviscosity parameter
nu_final = 3 * nu0 / 4;        % final hyperviscosity parameter

% Initialize continuation vectors / matrices
nu_vec = linspace(nu0, nu_final, m + 1);    % vector of hyperviscosity values
v_vec = zeros(round(N / 3), m + 1);         % matrix of converged states
T_vec = zeros(1, m + 1);                    % vector of converged periods
dt_vec = 0.1 * nu_vec;                      % vector of time-steps for time integration
lambda_vec = zeros(1, m + 1);               % vector of max. Floquet multiplier
flag_vec = zeros(1, m + 1);                 % vector of convergence flags

%% Converge and plot the limit cycle
% Initial condition at nu0
[x,~] = domain(N);
u_init = sin(2*pi*x) / sqrt(nu0);
v_init = field2vector(u_init,N);

% Transient integration to reach the limit cycle
T_trans = 1000 * nu0;
T_int = 250  * nu0;
dt = dt_vec(1);
[v_trans,~] = KSE_integrate(v_init,T_trans,dt,0,nu0,N);

% Store a short trajectory on the limit cycle
[v_traj,t_traj] = KSE_integrate(v_trans,T_int,dt,dt,nu0,N);

% Use the last state as initial guess for a point on the periodic orbit
vg0 = v_traj(:,end);

% Estimate the period from consecutive maxima of the energy
[E_traj,~,~] = projection(v_traj, N, nu0);
E_traj = E_traj(:);
t_traj = t_traj(:);

imax = find( ...
    E_traj(2:end-1) > E_traj(1:end-2) & ...
    E_traj(2:end-1) > E_traj(3:end) ...
) + 1;

Tg0 = t_traj(imax(end)) - t_traj(imax(end-1));

% Converge to the periodic orbit
[v0,T0,flag] = search4PO(vg0,Tg0,dt,nu0,N);

% Save converged solution
v_vec(:,1) = v0;
T_vec(1) = T0;
flag_vec(1) = flag;

% Verify by integrating the converged state for one period
[v_orbit,t_orbit] = KSE_integrate(v0,T0,dt,dt,nu0,N);

% Space-time plot
u_orbit = zeros(N,length(t_orbit));

for j = 1:length(t_orbit)
    u_orbit(:,j) = vector2field(v_orbit(:,j),N);
end

plot2C = figure("Name", "Plot2C");

contourf(x,t_orbit,u_orbit.',40,'LineColor','none');

xlabel('$x\;\mathrm{[m]}$','Interpreter','latex', ...
    'FontSize', params.plot.fontSize);

ylabel('$t\;\mathrm{[s]}$','Interpreter','latex', ...
    'FontSize', params.plot.fontSize);

colorbar;

% Save the plot automatically
if isgraphics(plot2C)
    exportgraphics(plot2C, fullfile(params.plot.folder, "Plot2C.png"), "Resolution", 300);
else
    disp('Skipped exporting plot2C (figure was closed).');
end

% 3D projection in (E,P,D)
[E_orbit,P_orbit,D_orbit] = projection(v_orbit,N,nu0);

plot2D = figure("Name", "Plot2D");
plot3(E_orbit, P_orbit, D_orbit,'LineWidth', params.plot.lw);
grid on;

xlabel('$E$','Interpreter','latex', ...
    'FontSize', params.plot.fontSize);

ylabel('$P$','Interpreter','latex', ...
    'FontSize', params.plot.fontSize);

zlabel('$D$','Interpreter','latex', ...
    'FontSize', params.plot.fontSize);

% Save the plot automatically
if isgraphics(plot2D)
    exportgraphics(plot2D, fullfile(params.plot.folder, "Plot2D.png"), "Resolution", 300);
else
    disp('Skipped exporting plot2D (figure was closed).');
end

%% Continuation
% First continuation guess from the converged limit cycle at nu0
i = 2;

nu = nu_vec(i);
dt = dt_vec(i);

vg = v_vec(:,i-1);
Tg = T_vec(i-1);

% Converge the periodic orbit at the new value of nu
[v_new,T_new,flag] = search4PO(vg,Tg,dt,nu,N);

% Save converged solution
v_vec(:,i) = v_new;
T_vec(i) = T_new;
flag_vec(i) = flag;

%% Subsequent continuation steps
i = 3;

while i <= m + 1
    nu = nu_vec(i);
    dt = dt_vec(i);

    % Linear extrapolation from the two previous converged solutions,
    % assuming uniform spacing in nu
    vg = 2 * v_vec(:,i-1) - v_vec(:,i-2);
    Tg = 2 * T_vec(i-1)   - T_vec(i-2);

    % Converge the periodic orbit at the new value of nu
    [v_new,T_new,flag] = search4PO(vg,Tg,dt,nu,N);

    % Save converged solution
    v_vec(:,i) = v_new;
    T_vec(i) = T_new;
    flag_vec(i) = flag;

    i = i + 1;
end

% Plot period as a function of nu
plot2E = figure('Name', 'Plot2E');
hold on; grid on;

% Identify successful and unsuccessful convergence
conv_idx = flag_vec > 0;
fail_idx = flag_vec <= 0;

% Plot continuation curve
p1 = plot(nu_vec, T_vec,'-','LineWidth', params.plot.lw);
handles = p1; 
labels = {'T curve'};

% Plot converged points
if any(conv_idx)
    p2 = plot(nu_vec(conv_idx), T_vec(conv_idx), 'o', ...
        'MarkerSize', 6, 'LineWidth', params.plot.lw);

    handles(end+1) = p2;
    labels{end+1} = 'converged';
end

% Plot non-converged points
if any(fail_idx)
    p3 = plot(nu_vec(fail_idx), T_vec(fail_idx), 'x', ...
        'MarkerSize', 8, 'LineWidth', params.plot.lw);

    handles(end+1) = p3;
    labels{end+1} = 'not converged';
end

xlabel('$\nu\;\mathrm{[m^2/s]}$','Interpreter','latex', ...
    'FontSize', params.plot.fontSize);

ylabel('$T\;\mathrm{[s]}$','Interpreter','latex', ...
    'FontSize', params.plot.fontSize);

legend(handles, labels, 'Location','northwest','Interpreter','latex',...
    'FontSize', params.plot.legendFontSize);

% Save the plot automatically
if isgraphics(plot2E)
    exportgraphics(plot2E, fullfile(params.plot.folder, "Plot2E.png"), "Resolution", 300);
else
    disp('Skipped exporting plot2E (figure was closed).');
end

%% Floquet exponents
for i = 1:(m + 1)
    nu = nu_vec(i);
    dt = dt_vec(i);
    T  = T_vec(i);
    v0 = v_vec(:,i);

    J = Jacobian(v0, T, epsilon, dt, nu, N);
    lambda = eig(J);

    % Remove the trivial multiplier lambda = 1
    lambda_tol = 1e-6;
    idx_trivial = find(abs(lambda - 1) < lambda_tol);
    
    % If lambda = 1 is not exact because of numerical error,
    % remove the closest multiplier to 1
    if isempty(idx_trivial)
        [~,idx_trivial] = min(abs(lambda - 1));
    
    end

    lambda(idx_trivial) = [];

    % Store the largest non-trivial multiplier modulus
    lambda_vec(i) = max(abs(lambda));
end

% Plot the largest non-trivial Floquet multiplier as a function of nu
plot2F = figure('Name', 'Plot2F');
hold on; grid on;

p1_F = plot(nu_vec, lambda_vec, '-', 'LineWidth', params.plot.lw);
handles_F = p1_F; 
labels_F = {'$|\lambda^\ast|$ curve'};

if any(conv_idx)
    p2_F = plot(nu_vec(conv_idx), lambda_vec(conv_idx), 'o', ...
        'MarkerSize', 6, 'LineWidth', params.plot.lw);
    handles_F(end+1) = p2_F;
    labels_F{end+1} = 'converged';
end

if any(fail_idx)
    p3_F = plot(nu_vec(fail_idx), lambda_vec(fail_idx), 'x', ...
        'MarkerSize', 8, 'LineWidth', params.plot.lw);
    handles_F(end+1) = p3_F;
    labels_F{end+1} = 'not converged';
end

yline(1,'--k', '$|\lambda^\ast| = 1$', ...
     'Interpreter', 'latex', 'fontsize', params.plot.fontSize, ...
     'LineWidth', params.plot.lw, ...
     'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off', ...
     'LabelHorizontalAlignment', 'right');

xlabel('$\nu\;\mathrm{[m^2/s]}$','Interpreter','latex', ...
    'FontSize', params.plot.fontSize);

ylabel('$|\lambda^\ast|$','Interpreter','latex', ...
    'FontSize', params.plot.fontSize);

legend(handles_F, labels_F, 'Location','northwest','Interpreter','latex',...
    'FontSize', params.plot.legendFontSize);

% Save the plot automatically
if isgraphics(plot2F)
    exportgraphics(plot2F, fullfile(params.plot.folder, "Plot2F.png"), "Resolution", 300);
else
    disp('Skipped exporting plot2F (figure was closed).');
end

% Estimate the crossing of |lambda*| = 1
idx_cross = find((lambda_vec(1:end-1)-1).*(lambda_vec(2:end)-1) <= 0,1,'first');

if ~isempty(idx_cross)
    nu_cross = nu_vec(idx_cross) + ...
        (1 - lambda_vec(idx_cross)) * ...
        (nu_vec(idx_cross+1) - nu_vec(idx_cross)) / ...
        (lambda_vec(idx_cross+1) - lambda_vec(idx_cross));
end

%% --- Display results ---
fprintf('================ Part 2.3 results ================\n');

selected_idx = [1 6 11 16 21];
fprintf('\n      nu              T              |lambda*|              flag\n');
for j = selected_idx
    fprintf(' %e    %e    %e    %d\n', nu_vec(j), T_vec(j), lambda_vec(j), flag_vec(j));
end

if ~isempty(idx_cross)
    fprintf('\nEstimated crossing of |lambda*| = 1:\n');
    fprintf('nu_cross = %.6e\n', nu_cross);
end