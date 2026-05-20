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

% Parameters based on parts 2 and 4
params.lambda = 0.01124216;
params.Lc = 0.3573;
params.eta = 0.000185;    

[u, U] = load_data(params.fileName, Inf);
dx = abs(U) / params.samplingFreq;

lmin = max(1e-4, 2*dx);
lmax = 1000;
l_S2 = logspace(log10(lmin), log10(lmax), 200)';

%% --- Flatness and Hyper-flatness ---
S2 = structure_function(u, l_S2, params.samplingFreq, U, 2);
S4 = structure_function(u, l_S2, params.samplingFreq, U, 4);
S6 = structure_function(u, l_S2, params.samplingFreq, U, 6);

smooth_window = 7;

F_raw = S4./S2.^2;
H_raw = S6./S2.^3;
F_smooth = smoothdata(F_raw, 'movmean', smooth_window, 'omitnan');
H_smooth = smoothdata(H_raw, 'movmean', smooth_window, 'omitnan');

%% --- Calibrate dimensional parameter D of the models ---
% Fit the observed slopes of F(l) and H(l) in the inertial range
fit_range = (l_S2 > params.lambda) & (l_S2 < params.Lc) & ...
            isfinite(F_smooth) & F_smooth > 0 & ...
            isfinite(H_smooth) & H_smooth > 0;

pF = polyfit(log(l_S2(fit_range)), log(F_smooth(fit_range)), 1);
pH = polyfit(log(l_S2(fit_range)), log(H_smooth(fit_range)), 1);

slope_F = pF(1);
slope_H = pH(1);

% Beta-model: F ~ l^(-C), H ~ l^(-2C)
C_beta_opt = -(slope_F + 2*slope_H) / 5;
D_beta_opt = 3 - C_beta_opt;

% Bifractal model: F ~ l^(-C/3), H ~ l^(-C)
C_bifractal_opt = -(3*slope_F + 9*slope_H) / 10;
D_bifractal_opt = 3 - C_bifractal_opt;

% Use calibrated values
D_beta = D_beta_opt;
D_bifractal = D_bifractal_opt;

C_beta = 3 - D_beta;
C_bifractal = 3 - D_bifractal;

% Model exponents
% Beta-model: all orders follow beta-model scaling
params.beta_E_exp = -(5/3 + C_beta/3);
F_beta_exp = -C_beta;
H_beta_exp = -2*C_beta;

% Bifractal model: p <= 3 follows K41, p >= 3 follows beta-model
params.bifractal_E_exp = -5/3;
F_bifractal_exp = -C_bifractal/3;
H_bifractal_exp = -C_bifractal;

%% --- Plot Flatness ---
plot1H = figure("Name", "Plot1H");

% Checking validity
valid_F_raw = isfinite(l_S2) & isfinite(F_raw) & F_raw > 0;
valid_F_smooth = isfinite(l_S2) & isfinite(F_smooth) & F_smooth > 0;

% Plot raw data in grey
loglog(l_S2(valid_F_raw), F_raw(valid_F_raw), ...
    'Color', [0.8 0.8 0.8], ...
    'DisplayName', '$F(l)$ (raw)');

hold on; grid on;

% Plot smoothed data as the main curve
loglog(l_S2(valid_F_smooth), F_smooth(valid_F_smooth), ...
    'LineWidth', params.plot.lw, ...
    'DisplayName', '$F(l)$ (smoothed)');

% Averaging the last values of F to get and estimate for the limit
F_large_scale = mean(F_smooth(end-5:end), 'omitnan');

yline(F_large_scale, '--k', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.insetFontSize, ...
    'Label', '$F(l\rightarrow\infty)$', ...
    'LabelHorizontalAlignment', 'right', ...
    'DisplayName', '$F(l\rightarrow\infty)$');

xlabel('$l\;[m]$', 'Interpreter', 'latex', ...
    'fontsize', params.plot.fontSize);

ylabel('$F(l) = \frac{S_4(l)}{S_2(l)^2}$', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.fontSize);

legend('show', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.legendFontSize, ...
    'Location', 'NorthEast');

xline(params.eta, '--k', '$\eta$', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.insetFontSize, ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(params.lambda, '--r', '$\lambda$', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.insetFontSize, ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(params.Lc, '--b', '$L_c$', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.insetFontSize, ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

% Add intermittency predictions to Plot H
l_pred = logspace(log10(params.lambda), log10(params.Lc), 100);

l_ref = sqrt(l_pred(1) * l_pred(end));
[~, idx_ref] = min(abs(l_S2 - l_ref));
F_ref = F_smooth(idx_ref);

loglog(l_pred, F_ref * (l_pred/l_ref).^F_beta_exp, '-.', ...
    'LineWidth', params.plot.lw, ...
    'Color', [0, 0.5, 0], ...
    'DisplayName', sprintf('$\\beta$-model ($\\propto l^{%.3f}$)', F_beta_exp));

loglog(l_pred, F_ref * (l_pred/l_ref).^F_bifractal_exp, ':', ...
    'LineWidth', params.plot.lw, ...
    'Color', 'b', ...
    'DisplayName', sprintf('bifractal ($\\propto l^{%.3f}$)', F_bifractal_exp));

legend('show', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.legendFontSize, ...
    'Location', 'NorthEast');

% Save the plot automatically
if isgraphics(plot1H)
    exportgraphics(plot1H, fullfile(params.plot.folder, "plot1H.png"), "Resolution", 300);
else
    disp('Skipped exporting plot1H (figure was closed).');
end

%% --- Correction to energy spectrum ---
% If PART3 has been previously run, the plot will be reopened using
% the .fig matlab file

plot1C_fig  = fullfile(params.plot.folder, "plot1C.fig");
plot1C_data = fullfile(params.plot.folder, "plot1C_data.mat");

if ~exist(plot1C_fig, "file")
    error("plot1C.fig not found. Run PART3_energy_spectrum.m first\n");
end

if ~exist(plot1C_data, "file")
    error("plot1C_data.mat not found. Run PART3_energy_spectrum.m first\n");
end

load(plot1C_data, "k", "Ek_smooth");
figs = findall(0, "Type", "figure", "Name", "Plot1C");

if ~isempty(figs)
    plot1C = figs(1);
else
    plot1C = openfig(plot1C_fig, "visible");
    plot1C.Name = "Plot1C";
end

figure(plot1C);
ax_all = findall(plot1C, "Type", "axes");
ax = ax_all(1);
hold(ax, "on");

delete(findall(ax, "Tag", "beta_E_part8"));
delete(findall(ax, "Tag", "bifractal_E_part8"));

% Old hardcoded dimensionality parameters, now they are calibrated for each
% model
% params.beta_E_exp = -51/30;
% params.bifractal_E_exp = -53/30;

k_anchor_idx = floor(length(k)/2000); 
k_pred = [k(k_anchor_idx), k(floor(length(k)/10))];

E_beta = Ek_smooth(k_anchor_idx) * (k_pred / k_pred(1)).^params.beta_E_exp;
E_bifractal = Ek_smooth(k_anchor_idx) * (k_pred / k_pred(1)).^params.bifractal_E_exp;

loglog(ax, k_pred, E_beta, '-.', ...
    'LineWidth', params.plot.lw, ...
    'DisplayName', sprintf('$\\beta$-model ($\\propto k^{%.2f}$)', params.beta_E_exp), ...
    'Color', [0, 0.5, 0], ...
    'Tag', 'beta_E_part8');

loglog(ax, k_pred, E_bifractal, ':', ...
    'LineWidth', params.plot.lw, ...
    'DisplayName', sprintf('bifractal ($\\propto k^{%.2f}$)', params.bifractal_E_exp), ...
    'Color', 'b', ...
    'Tag', 'bifractal_E_part8');

legend(ax, 'show', ...
    'Interpreter', 'latex', ...
    'FontSize', params.plot.legendFontSize, ...
    'Location', 'SouthWest');

% Save the plot automatically
if isgraphics(plot1C)
    savefig(plot1C, plot1C_fig);
    exportgraphics(plot1C, fullfile(params.plot.folder, "plot1C.png"), ...
    "Resolution", 300);
else
    disp('Skipped exporting plot1C (figure was closed).');
end

%% --- Plot Hyper-flatness ---
plot1I = figure("Name", "Plot1I");

% Check variables
valid_H_raw = isfinite(l_S2) & isfinite(H_raw) & H_raw > 0;
valid_H_smooth = isfinite(l_S2) & isfinite(H_smooth) & H_smooth > 0;

% Plot raw data in grey
loglog(l_S2(valid_H_raw), H_raw(valid_H_raw), ...
    'Color', [0.8 0.8 0.8], ...
    'DisplayName', '$H(l)$ (raw)');

hold on; grid on;

% Plot smoothed data 
loglog(l_S2(valid_H_smooth), H_smooth(valid_H_smooth), ...
    'LineWidth', params.plot.lw, ...
    'DisplayName', '$H(l)$ (smoothed)');

% Add intermittency corrections to plot I

% Predictions are based on the calibrated dimensions D_beta and D_bifractal.
% beta-model:      H(l) ~ l^(-2*C_beta)
% bifractal model: H(l) ~ l^(-C_bifractal)

% Choose a representative range 
l_pred = logspace(log10(params.lambda), log10(params.Lc), 100);

% Anchor the prediction lines 
l_ref = sqrt(l_pred(1) * l_pred(end));
[~, idx_ref] = min(abs(l_S2 - l_ref));
H_ref = H_smooth(idx_ref);

loglog(l_pred, H_ref * (l_pred / l_ref).^H_beta_exp, '-.', ...
    'LineWidth', params.plot.lw, ...
    'Color', [0, 0.5, 0], ...
    'DisplayName', sprintf('$\\beta$-model ($\\propto l^{%.3f}$)', H_beta_exp));

loglog(l_pred, H_ref * (l_pred / l_ref).^H_bifractal_exp, ':', ...
    'LineWidth', params.plot.lw, ...
    'Color', 'b', ...
    'DisplayName', sprintf('bifractal ($\\propto l^{%.3f}$)', H_bifractal_exp));

% Averaging the last values of F to get and estimate for the limit
H_large_scale = mean(H_smooth(end-5:end), 'omitnan');

yline(H_large_scale, '--k', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.insetFontSize, ...
    'Label', '$H(l\rightarrow\infty)$', ...
    'LabelHorizontalAlignment', 'right', ...
    'DisplayName', '$H(l\rightarrow\infty)$');

xlabel('$l\;[\mathrm{m}]$', 'Interpreter', 'latex', 'fontsize', ...
    params.plot.fontSize);

ylabel('$H(l) = \frac{S_6(l)}{S_2(l)^3}$', 'Interpreter', 'latex', ...
    'fontsize', params.plot.fontSize);

legend('show', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.legendFontSize, ...
    'Location', 'NorthEast');

xline(params.eta, '--k', '$\eta$', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.insetFontSize, ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(params.lambda, '--r', '$\lambda$', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.insetFontSize, ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(params.Lc, '--b', '$L_c$', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.insetFontSize, ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

% Save the plot automatically
if isgraphics(plot1I)
    exportgraphics(plot1I, fullfile(params.plot.folder, "plot1I.png"), "Resolution", 300);
else
    disp('Skipped exporting plot1I (figure was closed).');
end


%% --- Display results ---
fprintf('================ Part 1.8 results ================\n');
fprintf('Large-scale flatness estimate F_inf: %e\n', F_large_scale);
fprintf('Gaussian flatness reference: 3\n');
fprintf('Large-scale hyper-flatness estimate H_inf: %e\n', H_large_scale);
fprintf('Gaussian hyper-flatness reference: 15\n');

fprintf('\n[Calibrated dimensions]\n');
fprintf('beta-model:      D = %.6f, C = %.6f\n', D_beta, C_beta);
fprintf('bifractal model: D = %.6f, C = %.6f\n', D_bifractal, C_bifractal);

fprintf('\n[Scalings for results table]\n');
fprintf('beta-model:      E(k) ~ k^(%.6f)\n', params.beta_E_exp);
fprintf('bifractal model: E(k) ~ k^(%.6f)\n', params.bifractal_E_exp);
fprintf('beta-model:      H(l) ~ l^(%.6f)\n', H_beta_exp);
fprintf('bifractal model: H(l) ~ l^(%.6f)\n', H_bifractal_exp);

fprintf('\n[Additional flatness scalings]\n');
fprintf('beta-model:      F(l) ~ l^(%.6f)\n', F_beta_exp);
fprintf('bifractal model: F(l) ~ l^(%.6f)\n', F_bifractal_exp);