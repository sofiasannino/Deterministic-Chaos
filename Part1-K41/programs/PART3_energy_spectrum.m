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
params.fileName = 'veldata.txt';
params.samplingFreq = 20000;
params.plot.fontSize = 20;
params.plot.legendFontSize = 18;
params.plot.insetFontSize = 16;
params.plot.lw = 2.1;
params.plot.lwsmall = 1.8;

[u, U] = load_data(params.fileName, Inf);
params.U = U;

params.N = length(u);
dl = params.U / params.samplingFreq;
l = (0:params.N-1)' * dl;
params.L = l(end);

%% --- Energy spectrum ---
% Use the analytical normalization prefactor
spectrumNorm = (2 * pi * params.N^2) / params.L;
Ek = spectral_energy_density(u) / spectrumNorm;
Ek = Ek(:);

% Define wave numbers k
% A spectral feature at k corresponds to the physical wavelength l = 2*pi/k
dk = 2 * pi / params.L;
k = (0:(params.N/2))' * dk;

% Define E(k), k > 0
Ek(2:end-1) = 2 * Ek(2:end-1); 
k = k(2:end);
Ek = Ek(2:end);

% Signal filtering
% We use a moving average to smooth the signal
Ek_smooth = smoothdata(Ek, 'movmean', 500);

%% Plot 1C
plot1C = figure("Name", "Plot1C");
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
params.Lint_e = 2 * pi / kint_e;

% Labels
xlabel('$k \;[\mathrm{m^{-1}}]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize)

ylabel('$E(k) \; [\mathrm{m^3/s^2}]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize)

% Plot correlation length scales Lci
xline(kint_e, '--r', '$k_{L,\mathrm{int}}$', ...
    'Interpreter', 'latex', ...
    'FontSize', params.plot.insetFontSize, ...
    'Color', 'red', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(keta_e, '--', '$k_{\eta_e}$', ...
    'Interpreter', 'latex', ...
    'FontSize', params.plot.insetFontSize, ...
    'Color', 'magenta', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'right', ...
    'HandleVisibility', 'off');

xline(kLc, '--g', '$k_{L_c}$', ...
    'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
    'Color', 'blue', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

legend('Interpreter', 'latex', ...
       'FontSize', params.plot.legendFontSize, ...
       'Location', 'southwest');

% Save the plot automatically
% Also save the plot as image and MATLAB figure to be reused in PARTS 4 and 8
if isgraphics(plot1C)
    exportgraphics(plot1C, fullfile(params.plot.folder, "plot1C.png"), "Resolution", 300);
    savefig(plot1C, fullfile(params.plot.folder, "plot1C.fig"));
else
    disp('Skipped exporting plot1C (figure was closed).');
end

% Save data to be used in PART8
save(fullfile(params.plot.folder, "plot1C_data.mat"), ...
    "k", "Ek", "Ek_smooth");

%% Parseval's theorem
% 1/2 <u^2> = int E(k) dk
LHS = 0.5 * mean(u.^2);
RHS = sum(Ek) * dk;
rel_error = abs(LHS - RHS) / LHS;

%% Additional plot to better estimate K41 scales
plot1Extra1= figure("Name", "Plot1Extra1");
loglog(k, k.^(5/3).*Ek_smooth, ...
    'LineWidth', params.plot.lw, ...
    'DisplayName', '$k^{5/3}E(k)$');

hold on; grid on;

xlabel('$k\;[\mathrm{m^{-1}}]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

ylabel('$k^{5/3}E(k)$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

title('Compensated spectrum', 'Interpreter', 'latex');

% Add K41 scaling predictions to the compensated spectrum
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

xline(kLc, '--b', '$k_{L_c}$', ...
    'Interpreter', 'latex', ...
    'FontSize', params.plot.insetFontSize, ...
    'Color', 'blue', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

legend('Interpreter', 'latex', ...
    'FontSize', params.plot.legendFontSize, ...
    'Location', 'northwest');

% Save the plot automatically
if isgraphics(plot1Extra1)
    exportgraphics(plot1Extra1, fullfile(params.plot.folder, "plot1Extra1.png"), "Resolution", 300);
else
    disp('Skipped exporting plot1Extra1 (figure was closed).');
end

%% --- Display results ---
fprintf('================ Part 1.3 results ================\n');
fprintf('Spectrum norm: %e\n', spectrumNorm);
fprintf('Analytical prefactor: L/(2 pi N^2)\n');
fprintf('Numerical prefactor: %e\n', 1/spectrumNorm);
fprintf('Eta_e: %e\n', params.eta_e);
fprintf('Lint_e %e\n', params.Lint_e);
fprintf('Relative error from Parsevals theorem: %e\n', rel_error);