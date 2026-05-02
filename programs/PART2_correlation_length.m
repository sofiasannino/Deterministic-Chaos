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


%% --- Correlation length ---
params.U = U;
dl = params.U / params.samplingFreq;

lmax_autocorr = 1000;
Cl = autocorrelation(u, dl, lmax_autocorr);


% Length vectors
l = 0:dl:(dl*(length(Cl) - 1));

% Find intersections
Lc_idx = find(Cl <= exp(-1), 1);
params.Lc = l(Lc_idx);
disp('Lc');
disp(params.Lc);


%% --- Plot B ---
lmax_plot = 1;
lmax_idx = find(l >= lmax_plot, 1);
if isempty(lmax_idx)
    lmax_idx = length(l);
end

plotB = figure();
hold on; grid on;

plot(l(1:lmax_idx), Cl(1:lmax_idx), 'LineWidth', params.plot.lw);

xlabel('$l\;[m]$', 'Interpreter', 'latex', 'FontSize', params.plot.fontSize);
ylabel('$C(l)$', 'Interpreter', 'latex', 'FontSize', params.plot.fontSize);

% Add Lc to plots
xline(params.Lc, '--r', '$L_c$', ...
    'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

%% --- Integral scale ---
zero_cross_idx = find(Cl <= 0, 1);
if isempty(zero_cross_idx)
    zero_cross_idx = length(Cl);
end

% Compute the integral length scale using numerical integration
params.Lint = trapz(l(1:zero_cross_idx), Cl(1:zero_cross_idx));
disp('Lint')
disp(params.Lint);
