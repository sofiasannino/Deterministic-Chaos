%% --- Clear workspace ---
clc;
clear;
addpath("../functions/");


%% --- Common setup ---
% Define any quantity you need for this part in here
params.fileName = 'veldata.txt';

% Sampling frequency is 20 kHz
params.samplingFreq = 20000;

params.plot.fontSize = 14;
params.plot.legendFontSize = 12;
params.plot.insetFontSize = 10;
params.plot.lw = 1.5;
params.plot.lwsmall = 1.2;

[u, U] = load_data(params.fileName, Inf);

%% --- Frozen flow hypothesis ---
% Create the time vector based on the number of samples and the sampling frequency
time_vector = (0:length(u)-1)' / params.samplingFreq;

params.U = U;

% Apply Taylor's frozen flow hypothesis: x = U * t
params.x = time_vector * params.U;
params.L = params.x(end);

%% --- Plot A ---
plotA = figure();
hold on; grid on;

% Plot total velocity u_tot(x) against x
plot(params.x, u + params.U, 'LineWidth', params.plot.lwsmall);

yline(params.U, '--r', 'Mean Velocity $U$', ...
     'Interpreter', 'latex', 'fontsize', params.plot.fontSize, ...
     'Color', 'white', ...
     'LineWidth', params.plot.lw, ...
     'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off', ...
     'LabelHorizontalAlignment', 'left');

xlabel('$x\;[m]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);
ylabel('$u_\mathrm{tot}(x)\;[m/s]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);


%% --- Turbulence intensity ---
params.I = sqrt(mean(u.^2)) / params.U;
disp(params.I)