%% --- Clear workspace ---
clc;
clear;
addpath("../functions/");


%% --- Common setup ---
% Define any quantity you need for this part in here
params.plot.fontSize = 14;
params.plot.legendFontSize = 12;
params.plot.insetFontSize = 10;
params.plot.lw = 1.5;
params.plot.lwsmall = 1.2;

[u, U] = load_data(params.fileName, Inf);


%% --- Correlation length ---
dl = 1;

lmax_autocorr = 1000;
Cl = autocorrelation(u, dl, lmax_autocorr);


% Length vectors
l = 0:dl:(dl*(length(Cl) - 1));

% Find intersections
Lc_idx = 1;
params.Lc = l(Lc_idx);

warning('Make sure to set the parameters properly properly')


%% --- Plot B ---
lmax_plot = 1;
lmax_idx = find(l >= lmax_plot, 1);
if isempty(lmax_idx)
    lmax_idx = length(l);
end

plotB = figure();
hold on; grid on;
% plot(..., 'LineWidth', params.plot.lw);

xlabel('$l\;[m]$', 'Interpreter', 'latex', 'FontSize', params.plot.fontSize);
ylabel('$C(l)$', 'Interpreter', 'latex', 'FontSize', params.plot.fontSize);

% Add Lc to plots
% xline(..., ...
%     'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
%     'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');


%% --- Integral scale ---
warning('Compute the integral length scale');
