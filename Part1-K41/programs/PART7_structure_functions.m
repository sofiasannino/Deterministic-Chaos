%% --- Move to this file's folder ---
cd(fileparts(mfilename('fullpath')));

%% --- Clear workspace ---
clc;
clear;
addpath("../functions/");

%% --- Output folder for plots ---
params.plot.folder = "../plots/part1";

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

% Added parameters based on parts 2 and 4
params.lambda = 0.01124216;
params.L_c = 0.3573;

[u, U] = load_data(params.fileName, Inf);

dx = abs(U) / params.samplingFreq;

lmin = max(1e-3, 2*dx);
lmax = 10;

l_S2 = logspace(log10(lmin), log10(lmax), 200);
l_S3 = l_S2;

% Scaling exponents
params.S2_exp  = 2/3;
params.S3_exp  = 1;
S2_exp_str = '2/3';
S3_exp_str = '1';

%% --- Compute S2 and S3 ---
S2 = structure_function(u, l_S2, params.samplingFreq, U, 2);
S3 = structure_function(u, l_S3, params.samplingFreq, U, 3);

%% --- Smooth data ---
smooth_window = 7;

S2_smooth = smoothdata(S2, 'movmean', smooth_window, 'omitnan');
S3_smooth = smoothdata(S3, 'movmean', smooth_window, 'omitnan');

%% --- Plot S2 ---
plotF = figure();

% Checking validity
valid_S2_raw = isfinite(S2) & S2 > 0;
valid_S2_smooth = isfinite(S2_smooth) & S2_smooth > 0;

% Plot the raw data in grey
loglog(l_S2(valid_S2_raw), S2(valid_S2_raw), ...
    'Color', [0.8 0.8 0.8], ...
    'DisplayName', '$S_2(l)$ (raw)');

hold on; grid on;

% Required loglog plot
loglog(l_S2(valid_S2_smooth), S2_smooth(valid_S2_smooth), ...
    'LineWidth', params.plot.lw, ...
    'DisplayName', '$S_2(l)$ (smoothed)');

hold on; grid on;

% Add K41 scaling prediction.
% K41 predicts S2(l) ~ l^(2/3) in the inertial range.
S2_anchor_idx = floor(length(l_S2)/2);
S2_k41 = S2_smooth(S2_anchor_idx) * ...
    (l_S2 / l_S2(S2_anchor_idx)).^params.S2_exp;

loglog(l_S2, S2_k41, 'k--', ...
    'LineWidth', params.plot.lwsmall, ...
    'DisplayName', ['K41 ($l^{', S2_exp_str, '}$)']);

xlabel('$l\;[m]$', 'Interpreter', 'latex', ...
    'fontsize', params.plot.fontSize);

ylabel('$S_2(l)\;[m^2/s^2]$', 'Interpreter', 'latex', ...
    'fontsize', params.plot.fontSize);

legend('show', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.legendFontSize, ...
    'Location','NorthWest');

xline(params.lambda, '--r', '$\lambda$',...
    'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

xline(params.L_c, '--b', '$L_c$', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.insetFontSize, ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

%% --- Plot S3 ---
plotG = figure();

% S3 can be negative because of the four-fifth law.
% For a log-log plot we need its magnitude.
S3_magnitude = abs(S3);
S3_magnitude_smooth = abs(S3_smooth);

valid_S3_raw = isfinite(S3_magnitude) & S3_magnitude > 0;
valid_S3_smooth = isfinite(S3_magnitude_smooth) & S3_magnitude_smooth > 0;

% Raw data in grey
loglog(l_S3(valid_S3_raw), S3_magnitude(valid_S3_raw), ...
    'Color', [0.8 0.8 0.8], ...
    'DisplayName', '$|S_3(l)|$ (raw)');

hold on; grid on;

loglog(l_S3(valid_S3_smooth), S3_magnitude_smooth(valid_S3_smooth), ...
    'LineWidth', params.plot.lw, ...
    'DisplayName', '$|S_3(l)|$ (smoothed)');

% Add K41 four-fifth law scaling prediction.
% K41 predicts |S3(l)| ~ l.
S3_anchor_idx = floor(length(l_S3)/2);
S3_k41 = S3_magnitude_smooth(S3_anchor_idx) * ...
    (l_S3 / l_S3(S3_anchor_idx)).^params.S3_exp;

loglog(l_S3, S3_k41, 'k--', ...
    'LineWidth', params.plot.lwsmall, ...
    'DisplayName', ['K41 ($l^{', S3_exp_str, '}$)']);

xlabel('$l\;[m]$', 'Interpreter', 'latex', ...
    'fontsize', params.plot.fontSize);

ylabel('$|S_3(l)|\;[m^3/s^3]$', 'Interpreter', 'latex', ...
    'fontsize', params.plot.fontSize);

legend('show', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.legendFontSize, ...
    'Location','NorthWest');

xline(params.lambda, '--r', '$\lambda$', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.insetFontSize, ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(params.L_c, '--b', '$L_c$', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.insetFontSize, ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

%% --- Estimate ε from S2 and S3 ---
C2 = 2.15;

% K41 estimates are meaningful in the inertial range so we restrict to the
% interval between the Taylor scale and the integral scale
lambda = params.lambda;
Lc = params.L_c;

lS2 = l_S2(:);
lS3 = l_S3(:);
S2sm = S2_smooth(:);
S3sm = S3_smooth(:);

inertial_S2 = (lS2 > lambda) & (lS2 < Lc) & ...
              isfinite(S2sm) & (S2sm > 0);

inertial_S3 = (lS3 > lambda) & (lS3 < Lc) & ...
              isfinite(S3sm);

% Estimate from S3:
% S3(l) = -(4/5) epsilon l
epsilon_S3_local = -(5/4) * S3sm ./ lS3;

valid_eps_S3 = inertial_S3 & isfinite(epsilon_S3_local) & ...
               (epsilon_S3_local > 0);

epsilon_S3 = median(epsilon_S3_local(valid_eps_S3), 'omitnan');

epsilon_S2_local = (S2sm / C2).^(3/2) ./ lS2;

valid_eps_S2 = inertial_S2 & isfinite(epsilon_S2_local) & ...
               (epsilon_S2_local > 0);

epsilon_S2 = median(epsilon_S2_local(valid_eps_S2), 'omitnan');

% Outer-scale estimate from Section 1.4, for comparison
epsilon_outer = 2.939219; 

fprintf('epsilon from S2:           %e\n', epsilon_S2);
fprintf('epsilon from S3:           %e\n', epsilon_S3);
fprintf('epsilon from outer scale:  %e\n', epsilon_outer);
fprintf('relative diff S2/outer:    %.2f%%\n', ...
        100 * abs(epsilon_S2 - epsilon_outer) / epsilon_outer);
fprintf('relative diff S3/outer:    %.2f%%\n', ...
        100 * abs(epsilon_S3 - epsilon_outer) / epsilon_outer);

%% --- Add result on plots F and G ---
S2_epsilon_fit = C2 * (epsilon_S2 .* lS2).^(2/3);

% Plot the magnitude of the four-fifth law, since Plot G is loglog of |S3|
S3_epsilon_fit = (4/5) * epsilon_S3 .* lS3;

figure(plotF);
hold on
loglog(lS2, S2_epsilon_fit,    'LineWidth', params.plot.lw, 'DisplayName', ...
    '$y = C_2 \left( \langle \epsilon \rangle l \right)^{2/3}$');

figure(plotG);
hold on;

loglog(lS3, S3_epsilon_fit,...
    'LineWidth', params.plot.lw, 'DisplayName', ...
    '$y = \frac{4}{5} \langle \epsilon \rangle l $');

% Save the plot automatically
exportgraphics(plotF, fullfile(params.plot.folder, "plotF.png"), ...
    "Resolution", 300);
exportgraphics(plotG, fullfile(params.plot.folder, "plotG.png"), ...
    "Resolution", 300);
