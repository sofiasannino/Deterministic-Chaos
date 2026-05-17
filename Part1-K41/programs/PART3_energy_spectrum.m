%% --- Move to this file's folder ---
cd(fileparts(mfilename('fullpath')));

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
params.N = N;
dx = abs(params.U) / params.samplingFreq;
L = N * dx;
params.L = L;

%% --- Energy spectrum ---
% Use the analytically derived normalization prefactor
spectrumNorm = (2 * pi * params.N^2) / params.L;
Ek = spectral_energy_density(u) / spectrumNorm;
Ek = Ek(:);
fprintf('Spectrum norm: %e\n', spectrumNorm);
fprintf('Value for the table, analytical prefactor: L/(2 pi N^2)\n');
fprintf('Value for the table, numerical prefactor: %e\n', 1/spectrumNorm);

% Define wave numbers k
% These are angular wavenumbers, coherent with the Fourier convention
dk = 2*pi / L;
k = (0:(N/2))' * dk;

% exp(-ikx) used in the assignment. Therefore, a spectral feature at k
% corresponds to the physical wavelength l = 2*pi/k. The relation k ~ 1/l
% is only an order-of-magnitude scaling convention

% Define E(k), k > 0
Ek(2:end-1) = 2 * Ek(2:end-1); 
k = k(2:end);
Ek = Ek(2:end);

% Signal filtering
% We use a moving average to smooth the signal
Ek_smooth = smoothdata(Ek, 'movmean', 500);

%% Plot C
plotC = figure("Name", "PlotC");
loglog(k, Ek, 'Color', [0.8 0.8 0.8], 'DisplayName', '$E(k)$ (raw)');
hold on; grid on;

loglog(k, Ek_smooth, 'DisplayName', '$E(k)$ (smoothed)', ...
     'LineWidth', params.plot.lw);

% Add K41 scaling predictions
% K41 predicts E(k) scales as k^{-5/3} in the inertial subrange
k_anchor_idx = floor(length(k)/2000); 
k41_k = [k(k_anchor_idx), k(floor(length(k)/10))];
k41_E = Ek_smooth(k_anchor_idx) * (k41_k / k41_k(1)).^(-5/3);
loglog(k41_k, k41_E, 'k--', 'LineWidth', params.plot.lw, 'DisplayName', 'K41 ($k^{-5/3}$)');

% Graphical estimates for integral and Kolmogorov length scales
keta_e = 615;
kint_e = 5;
kLc = 2 * pi / 0.3573;

params.eta_e = 2 * pi / keta_e;
fprintf('Eta_e: %e\n', params.eta_e);

params.Lint_e = 2 * pi / kint_e;
fprintf('Lint_e %e\n', params.Lint_e);

% Labels
xlabel('$k \;[m^{-1}]$', 'Interpreter', 'latex')
ylabel('$E(k) \; [m^3/s^2]$', 'Interpreter', 'latex')


% Plot correlation length scales Lci
xline(kint_e, '--r', '$k_{L,\mathrm{int}}$', ...
    'Interpreter', 'latex', ...
    'FontSize', params.plot.insetFontSize, ...
    'Color', 'red', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'DisplayName', '$k_{L_{\mathrm{int}}}$');

xline(keta_e, '--', '$k_{\eta_e}$', ...
    'Interpreter', 'latex', ...
    'FontSize', params.plot.insetFontSize, ...
    'Color', 'magenta', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'DisplayName', '$k_{\eta,e}$');

xline(kLc, '--g', '$k_{L_c}$', ...
    'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
    'Color', 'blue', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left','HandleVisibility', 'off');

legend('Interpreter', 'latex', ...
       'FontSize', params.plot.legendFontSize, ...
       'Location', 'best');

% Save Plot C as image and MATLAB figure to be reused later in PART 4 and 8
plotC_png = fullfile(params.plot.folder, "plotC.png");
plotC_fig = fullfile(params.plot.folder, "plotC.fig");

exportgraphics(plotC, plotC_png, "Resolution", 300);
savefig(plotC, plotC_fig);

% Save data to be used in PART8
save(fullfile(params.plot.folder, "plotC_data.mat"), ...
    "k", "Ek", "Ek_smooth");

%% Parsival's theorem
% 1/2 <u^2> = int E(k) dk
LHS = 0.5 * mean(u.^2);
RHS = sum(Ek) * dk;
rel_error = abs(LHS - RHS) / LHS;
fprintf('Relative Error: %e\n', rel_error);

%% Additional plot to better estimate K41 scales
plotC_comp = figure("Name", "Compensated spectrum");

loglog(k, k.^(5/3).*Ek_smooth, ...
    'LineWidth', params.plot.lw, ...
    'DisplayName', '$k^{5/3}E(k)$');

hold on; grid on;

xlabel('$k\;[m^{-1}]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

ylabel('$k^{5/3}E(k)$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

title('Compensated spectrum', 'Interpreter', 'latex');

% --- K41 prediction on compensated spectrum ---
% If E(k) ~ k^(-5/3), then k^(5/3)E(k) is approximately constant
k41_comp_idx_1 = floor(length(k)/2000);
k41_comp_idx_2 = floor(length(k)/10);

k41_comp_k = [k(k41_comp_idx_1), k(k41_comp_idx_2)];
k41_comp_level = k(k41_comp_idx_1)^(5/3) * Ek_smooth(k41_comp_idx_1);

loglog(k41_comp_k, ...
       k41_comp_level * [1 1], ...
       'k--', ...
       'LineWidth', params.plot.lw, ...
       'DisplayName', 'K41 plateau');

% --- Characteristic spectral scales ---

xline(kint_e, '--r', '$k_{\mathrm{int},e}$', ...
    'Interpreter', 'latex', ...
    'FontSize', params.plot.insetFontSize, ...
    'Color', 'red', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(keta_e, '--g', '$k_{\eta,e}$', ...
    'Interpreter', 'latex', ...
    'FontSize', params.plot.insetFontSize, ...
    'Color', 'magenta', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(kLc, '--b', '$k_{L_c}=2\pi/L_c$', ...
    'Interpreter', 'latex', ...
    'FontSize', params.plot.insetFontSize, ...
    'Color', 'blue', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

legend('Interpreter', 'latex', ...
    'FontSize', params.plot.legendFontSize, ...
    'Location', 'northwest');

%% --- Save compensated spectrum ---
exportgraphics(plotC_comp, ...
    fullfile(params.plot.folder, "plotC_compensated.png"), ...
    "Resolution", 300);


