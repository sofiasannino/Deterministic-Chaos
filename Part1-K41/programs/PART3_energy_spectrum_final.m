%% --- Clear workspace ---
clc;
clear;
addpath("../functions/");


%% --- Common setup ---
% Define any quantity you need for this part in here
params.fileName = 'veldata.txt';
params.samplingFreq = 20000;
params.plot.fontSize = 14;
params.plot.legendFontSize = 12;
params.plot.insetFontSize = 10;
params.plot.lw = 1.5;
params.plot.lwsmall = 1.2;

[u, U] = load_data(params.fileName, Inf);

params.U = U;
dl = params.U / params.samplingFreq;
N = length(u);
l = (0:N-1)' * dl;
L = l(N);
params.N = N;
params.L = L;
Lc = 3.57303e-01;

%% --- Energy spectrum ---
% Use the analytically derived normalization prefactor
spectrumNorm = (2 * pi * params.N^2) / params.L;
Ek = spectral_energy_density(u./sqrt(spectrumNorm));
fprintf('Prefactor: %e\n', 1./spectrumNorm);

% Define wave numbers k
dk = 2 * pi / params.L;
k = (0:(N/2))' * dk;

% Define E(k), k > 0
Ek(2:end-1) = 2 * Ek(2:end-1); 
k = k(2:end);
Ek = Ek(2:end);

% Signal filtering
% We use a moving average to smooth the signal
Ek_smooth = smoothdata(Ek, 'movmean', 500);

%% Plot C
plotC = figure();
loglog(k, Ek, 'Color', [0.8 0.8 0.8], 'DisplayName', '$E(k)$ (raw)');
hold on; grid on;

loglog(k, Ek_smooth, 'DisplayName', '$E(k)$ (smoothed)', ...
     'LineWidth', params.plot.lw);

% Add K41 scaling predictions
% K41 predicts E(k) scales as k^{-5/3} in the inertial subrange
k_anchor_idx = floor(length(k)/3000); 
k41_k = [k(k_anchor_idx), k(floor(length(k)/10))];
k41_E = Ek_smooth(k_anchor_idx) * (k41_k / k41_k(1)).^(-5/3);
loglog(k41_k, k41_E, 'k--', 'LineWidth', params.plot.lw, 'DisplayName', 'K41 ($k^{-5/3}$)');

% Graphical estimates for integral and Kolmogorov length scales
keta_e = 615;
kint_e = 2;
kLc = 2 * pi / Lc;

params.eta_e = 2 * pi / keta_e;
fprintf('Eta_e: %e\n', params.eta_e);

params.Lint_e = 2 * pi / kint_e;
fprintf('Lint_e %e\n', params.Lint_e);

% Labels
xlabel('$k \;[m^{-1}]$', 'Interpreter', 'latex')
ylabel('$E(k) \; [m^3/s^2]$', 'Interpreter', 'latex')
legend('Interpreter', 'latex', 'FontSize', params.plot.legendFontSize);

% Plot correlation length scales Lci
xline(kint_e, '--r', '$k_{int}$', ...
    'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
    'Color', 'red', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left','HandleVisibility', 'off');
    
xline(keta_e, '--g', '$k_{\eta}$', ...
    'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
    'Color', 'red', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left','HandleVisibility', 'off');

xline(kLc, '--g', '$k_{L_c}$', ...
    'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
    'Color', 'blue', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left','HandleVisibility', 'off');

%% Parsival's theorem
% 1/2 <u^2> = int E(k) dk
LHS = 0.5 * mean(u.^2);
RHS = trapz(k, Ek);
rel_error = abs(LHS - RHS) / LHS;
fprintf('Relative Error: %e\n', rel_error);