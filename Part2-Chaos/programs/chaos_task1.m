close all; clear; clc

%% inputs
addpath('../functions/')       % folder containing functions
N = 128;                       % spatial resolution

nu_values = [1e-2, 1e-3, 8 * 1e-4, 6.5 * 1e-4, 2 * 1e-4];
nu = nu_values(5);             % hyperviscosity parameter

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
r = zeros(size(v1000));        % the perturbation vector memory allocation
for k = 1:length(r)
    eta = rand();              % random number in the interval (0,1)
    r(k) = epsilon * (2*eta - 1) * v1000(k);
end

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
figure('Position', [100, 100, 1200, 800], 'Name', sprintf('KSE Dynamics for nu = %g', nu));

% Panel 1: Space-time contour of u1(x,t)
subplot(2, 2, 1);
contourf(x, tt1/nu, u1', 50, 'LineStyle', 'none');
colorbar;
xlabel('x'); ylabel('t / \nu');
title('Panel 1: Unperturbed u_1(x,t)');

% Panel 2: Space-time contour of u2(x,t)
subplot(2, 2, 2);
contourf(x, tt2/nu, u2', 50, 'LineStyle', 'none');
colorbar;
xlabel('x'); ylabel('t / \nu');
title('Panel 2: Perturbed u_2(x,t)');

% Panel 3: Absolute difference |u1(x,t) - u2(x,t)|
subplot(2, 2, 3);
abs_diff = abs(u1 - u2);
contourf(x, tt1/nu, abs_diff', 50, 'LineStyle', 'none');
colorbar;
xlabel('x'); ylabel('t / \nu');
title('Panel 3: Difference |u_1 - u_2|');

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

xlabel('Energy Production Rate (P)'); 
ylabel('Energy Dissipation Rate (D)');
title('Panel 4: 2D Phase Portrait');
legend('Location', 'best');

%% Extra: compute divergence timescale
error_norm = vecnorm(vv1 - vv2);

% Plot the error
figure('Name', sprintf('Error for nu = %g', nu));
semilogy(tt1/nu, error_norm, 'k', 'LineWidth', 1.5);
hold on; grid on;

xlabel('Time (t / \nu)');
ylabel('|| v_1(t) - v_2(t) ||');
title('Divergence of perturbed trajectory');

mean_state_size = mean(vecnorm(vv1));
threshold = 0.1 * mean_state_size;
threshold_label = '10% Mean Attractor Size';

% Alternative
% threshold = 100 * error_norm(1);
% threshold_label = '100x Initial Error';

yline(threshold, 'r--', threshold_label, 'LabelHorizontalAlignment', 'left');
dev_idx = find(error_norm >= threshold, 1);

if ~isempty(dev_idx)
    t_deviation = tt1(dev_idx) / nu;
    plot(t_deviation, error_norm(dev_idx), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    fprintf('For nu = %g, significant deviation occurs at t/nu approx: %.1f\n', nu, t_deviation);
else
    fprintf('For nu = %g, no significant deviation occurred (system is stable).\n', nu);
end