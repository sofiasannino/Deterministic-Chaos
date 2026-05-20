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

[u, U] = load_data(params.fileName, Inf);
dx = abs(U) / params.samplingFreq;

lmin = max(1e-3, 2*dx);
lmax = 10;
l_S2 = logspace(log10(lmin), log10(lmax), 200);
l_S3 = l_S2;

% Scaling exponents
params.S2_exp  = 2/3;
params.S3_exp  = 1;

%% --- Compute S2 and S3 ---
S2 = structure_function(u, l_S2, params.samplingFreq, U, 2);
S3 = structure_function(u, l_S3, params.samplingFreq, U, 3);

% Smoothing
smooth_window = 7;
S2_smooth = smoothdata(S2, 'movmean', smooth_window, 'omitnan');
S3_smooth = smoothdata(S3, 'movmean', smooth_window, 'omitnan');

%% --- Plot S2 ---
plot1F = figure("Name", "Plot1F");

u_rms = sqrt(mean(u.^2));
nu = 1.51e-5;
params.epsilon = 0.5 * (u_rms^3) / params.Lc;
params.eta = (nu^3 / params.epsilon)^(1/4);

% Checking validity
valid_S2_raw = isfinite(S2) & S2 > 0;
valid_S2_smooth = isfinite(S2_smooth) & S2_smooth > 0;

% Plot the raw data in grey
loglog(l_S2(valid_S2_raw), S2(valid_S2_raw), ...
    'Color', [0.8 0.8 0.8], ...
    'DisplayName', '$S_2(l)$ (raw)');

hold on; grid on;

loglog(l_S2(valid_S2_smooth), S2_smooth(valid_S2_smooth), ...
    'LineWidth', params.plot.lw, ...
    'DisplayName', '$S_2(l)$ (smoothed)');

xlabel('$l\;[\mathrm{m}]$', 'Interpreter', 'latex', ...
    'fontsize', params.plot.fontSize);

ylabel('$S_2(l)\;[\mathrm{m^2/s^2}]$', 'Interpreter', 'latex', ...
    'fontsize', params.plot.fontSize);

legend('show', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.legendFontSize, ...
    'Location','NorthWest');

xline(params.lambda, '--r', '$\lambda$',...
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

xline(params.eta, '--k', '$\eta$', ... 
    'Interpreter', 'latex', ... 
    'fontsize', params.plot.insetFontSize, ... 
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ... 
    'HandleVisibility', 'off');

% Pope's dissipation-range estimate
% eta * 8 < l < eta * 60
l_eta_min = params.eta * 8;
l_eta_max = params.eta * 60;

p_eta_range = xregion(l_eta_min, l_eta_max, ...
    'FaceColor', [1.00 0.78 0.25], ...
    'FaceAlpha', 0.45, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');

uistack(p_eta_range, 'bottom');

xline(l_eta_max, '--', '$\eta$', ...
'Interpreter', 'latex', ...
'FontSize', params.plot.insetFontSize, ...
'Color', '#FF8C00', ...
'LabelVerticalAlignment', 'top', ...
'LabelHorizontalAlignment', 'left', ...
'DisplayName', '\eta$');

%% --- Plot S3 ---
plot1G = figure("Name", "Plot1G");

% S3 can be negative, we plot its magnitude
S3_magnitude = abs(S3);
S3_magnitude_smooth = abs(S3_smooth);

% Checking validity
valid_S3_raw = isfinite(S3_magnitude) & S3_magnitude > 0;
valid_S3_smooth = isfinite(S3_magnitude_smooth) & S3_magnitude_smooth > 0;

% Plot raw data in grey
loglog(l_S3(valid_S3_raw), S3_magnitude(valid_S3_raw), ...
    'Color', [0.8 0.8 0.8], ...
    'DisplayName', '$|S_3(l)|$ (raw)');

hold on; grid on;

loglog(l_S3(valid_S3_smooth), S3_magnitude_smooth(valid_S3_smooth), ...
    'LineWidth', params.plot.lw, ...
    'DisplayName', '$|S_3(l)|$ (smoothed)');

xlabel('$l\;[\mathrm{m}]$', 'Interpreter', 'latex', ...
    'fontsize', params.plot.fontSize);

ylabel('$|S_3(l)|\;[\mathrm{m^3/s^3}]$', 'Interpreter', 'latex', ...
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

xline(params.Lc, '--b', '$L_c$', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.insetFontSize, ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(params.eta, '--k', '$\eta$', ... 
    'Interpreter', 'latex', ... 
    'fontsize', params.plot.insetFontSize, ... 
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ... 
    'HandleVisibility', 'off');

% Pope's dissipation-range estimate
% eta * 8 < l < eta * 60
l_eta_min = params.eta * 8;
l_eta_max = params.eta * 60;

p_eta_range = xregion(l_eta_min, l_eta_max, ...
    'FaceColor', [1.00 0.78 0.25], ...
    'FaceAlpha', 0.45, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off');

uistack(p_eta_range, 'bottom');

xline(l_eta_max, '--', '$\eta$', ...
'Interpreter', 'latex', ...
'FontSize', params.plot.insetFontSize, ...
'Color', '#FF8C00', ...
'LabelVerticalAlignment', 'top', ...
'LabelHorizontalAlignment', 'left', ...
'DisplayName', '\eta$');

%% --- Estimate epsilon from S2 and S3 ---
% Sreenivasan's prefactor of the second order structure function
C2 = 2.15;

% K41 estimates are meaningful in the inertial range so we restrict to the
% interval between the Taylor scale and the integral scale
lS2 = l_S2(:);
lS3 = l_S3(:);
S2sm = S2_smooth(:);
S3sm = S3_smooth(:);

inertial_S2 = (lS2 > params.lambda) & (lS2 < params.Lc) & ...
              isfinite(S2sm) & (S2sm > 0);

inertial_S3 = (lS3 > params.lambda) & (lS3 < params.Lc) & ...
              isfinite(S3sm);

% S3(l) = -(4/5) epsilon l
epsilon_S3_local = -(5/4) * S3sm ./ lS3;

valid_eps_S3 = inertial_S3 & isfinite(epsilon_S3_local) & ...
               (epsilon_S3_local > 0);

epsilon_S3 = median(epsilon_S3_local(valid_eps_S3), 'omitnan');

% S2(l) = C2 (epsilon l)^(2/3)
epsilon_S2_local = (S2sm / C2).^(3/2) ./ lS2;

valid_eps_S2 = inertial_S2 & isfinite(epsilon_S2_local) & ...
               (epsilon_S2_local > 0);

epsilon_S2 = median(epsilon_S2_local(valid_eps_S2), 'omitnan');

% Outer-scale estimate from Section 1.4, for comparison
epsilon_outer = 2.939219; 

%% --- Add result on plots F and G ---
S2_epsilon_fit = C2 * (epsilon_S2 .* lS2).^(2/3);

% Plot the magnitude of the four-fifth law, since Plot G is loglog of |S3|
S3_epsilon_fit = (4/5) * epsilon_S3 .* lS3;

figure(plot1F);
hold on;

loglog(lS2, S2_epsilon_fit,    'LineWidth', params.plot.lw, 'DisplayName', ...
    '$C_2 \left( \langle \epsilon \rangle l \right)^{2/3}$');

yl_F = ylim;
text(sqrt(l_eta_min*l_eta_max), yl_F(1)*(yl_F(2)/yl_F(1))^0.08, ...
    '$\eta\cdot8 < l < \eta\cdot60$', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.insetFontSize, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom');

figure(plot1G);
hold on;

loglog(lS3, S3_epsilon_fit,...
    'LineWidth', params.plot.lw, 'DisplayName', ...
    '$\frac{4}{5} \langle \epsilon \rangle l $');

yl_G = ylim;
text(sqrt(l_eta_min*l_eta_max), yl_G(1)*(yl_G(2)/yl_G(1))^0.08, ...
    '$\eta\cdot8 < l < \eta\cdot60$', ...
    'Interpreter', 'latex', ...
    'fontsize', params.plot.insetFontSize, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom');

% Save the plots automatically
if isgraphics(plot1F)
    exportgraphics(plot1F, fullfile(params.plot.folder, "plot1F.png"), "Resolution", 300);
else
    disp('Skipped exporting plot1F (figure was closed).');
end

if isgraphics(plot1G)
    exportgraphics(plot1G, fullfile(params.plot.folder, "plot1G.png"), "Resolution", 300);
else
    disp('Skipped exporting plot1G (figure was closed).');
end

%% --- Display results ---
fprintf('================ Part 1.7 results ================\n');
fprintf('epsilon from S2: %e\n', epsilon_S2);
fprintf('epsilon from S3: %e\n', epsilon_S3);
fprintf('epsilon from outer scale: %e\n', epsilon_outer);
fprintf('relative diff S2/outer: %f%%\n', ...
        100 * abs(epsilon_S2 - epsilon_outer) / epsilon_outer);
fprintf('relative diff S3/outer: %f%%\n', ...
        100 * abs(epsilon_S3 - epsilon_outer) / epsilon_outer);