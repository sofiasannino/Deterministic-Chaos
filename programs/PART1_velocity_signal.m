%% --- Clear workspace ---
clc;
clear;
addpath("../functions/");


%% --- Common setup ---
% Define any quantity you need for this part in here
params.fileName = 'veldata.txt';
params.samplingFreq = 1;
params.plot.fontSize = 14;
params.plot.legendFontSize = 12;
params.plot.insetFontSize = 10;
params.plot.lw = 1.5;
params.plot.lwsmall = 1.2;

[u, U] = load_data(params.fileName, Inf);


%% --- Frozen flow hypothesis ---
params.U = 1;
params.x = 1;
params.L = 1;
warning('Make sure to properly fill the above parameters')


%% --- Plot A ---
plotA = figure();
hold on; grid on;
% plot(..., 'LineWidth', 1);

% yline(..., ...
%     'Interpreter', 'latex', 'fontsize', params.plot.fontSize, ...
%     'LineWidth', params.plot.lw, ...
%     'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off', ...
%     'LabelHorizontalAlignment', 'left');

xlabel('$x\;[m]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);
ylabel('$u_\mathrm{tot}(x)\;[m/s]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);


%% --- Turbulence intensity ---
params.I = 1;
warning('Make sure to fill the turbulence intensity value');