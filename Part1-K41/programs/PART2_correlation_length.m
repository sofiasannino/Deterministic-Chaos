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
params.plot.fontSize = 14;
params.plot.legendFontSize = 12;
params.plot.insetFontSize = 10;
params.plot.lw = 1.5;
params.plot.lwsmall = 1.2;

[u, U] = load_data(params.fileName, Inf);
params.U = U;

params.dx = - U./params.samplingFreq;
params.x = (0:length(u)-1)' * params.dx;
params.L = abs(params.x(end));

%% --- Correlation length ---
% Absolute displacement
dl = abs(params.dx);
lmax_autocorr = 1000;

% This function calculates correlation in accordance with our convention:
% a correlation of +l means we take all couples of velocities with index 
% distance freq * l / U and they are computed excluding the 'physical right
% tail' of the velocities, meaning that for an index q, velocities
% 0, 1, ..., q-1 (which physically are on the right of our domain, close to
% zero) won't be considered since we don't have data for velocities at
% physical distance +l from them
[Cl, normalization] = autocorrelation(u, dl, lmax_autocorr);

% Absolute length vector
l = 0:dl:(dl*(length(Cl) - 1));

% Find intersections
critic_value = 1/exp(1);
Lc_idx = find(Cl <= critic_value, 1, 'first');
params.Lc = l(Lc_idx);

%% --- Plot 1B ---
lmax_plot = 10;
lmax_idx = find(l >= lmax_plot, 1);
if isempty(lmax_idx)
    lmax_idx = length(l);
end

plot1B = figure("Name", "Plot1B");
hold on; grid on;
plot(l(1:lmax_idx), Cl(1:lmax_idx), 'LineWidth', params.plot.lw);

xlabel('$l\;[\mathrm{m}]$', 'Interpreter', 'latex', 'FontSize', params.plot.fontSize);
ylabel('$C(l)$', 'Interpreter', 'latex', 'FontSize', params.plot.fontSize);

% Add Lc to plots
xline(params.Lc, '--r', 'Correlation length Lc', ...
    'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% Save the plot automatically
if isgraphics(plot1B)
    exportgraphics(plot1B, fullfile(params.plot.folder, "plot1B.png"), "Resolution", 300);
else
    disp('Skipped exporting plot1B (figure was closed).');
end

%% --- Integral scale ---
% First zero crossing: Cl becomes negative
idx_zero1 = find(Cl <= 0, 1, 'first');

% Second zero crossing: Cl becomes positive again
idx_zero2 = idx_zero1 + find(Cl(idx_zero1+1:end) >= 0, 1, 'first');

% Minimum of Cl between first and second zero crossing
[~, idx_min_rel] = min(Cl(idx_zero1:idx_zero2));
idx_min = idx_zero1 + idx_min_rel - 1;

% Integrals
L_zero1 = trapz(l(1:idx_zero1), Cl(1:idx_zero1));
L_min   = trapz(l(1:idx_min),   Cl(1:idx_min));
L_zero2 = trapz(l(1:idx_zero2), Cl(1:idx_zero2));

% Errors
err_zero1 = abs(L_zero1 - params.Lc) / params.Lc;
err_min   = abs(L_min   - params.Lc) / params.Lc;
err_zero2 = abs(L_zero2 - params.Lc) / params.Lc;

%% --- Display results ---
fprintf('================ Part 1.2 results ================\n');
fprintf('Normalization: %e\n', normalization);
fprintf('Minimum value of Cl: %e \n', Cl(idx_min))
fprintf('Lc = %e\n', params.Lc)
fprintf('Lint = %e\n', L_zero1)
fprintf('Relative difference |Lint - Lc|/Lc = %e\n', err_zero1)

fprintf('\n================ Further experiments ================\n');
fprintf('[Stop at first zero crossing]\n')
fprintf('L_int: %e\n', L_zero1)
fprintf('Error: %e\n', err_zero1)

fprintf('\n[Stop at minimum]\n')
fprintf('L_int: %e\n', L_min)
fprintf('Error: %e\n', err_min)

fprintf('\n[Stop at second zero crossing]\n')
fprintf('L_int: %e\n', L_zero2)
fprintf('Error: %e\n', err_zero2)