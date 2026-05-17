%% --- Move to this file's folder ---
cd(fileparts(mfilename('fullpath')));

%% --- Clear workspace ---
clc;
clear;
addpath("../functions/");
addpath("../data/")

%% --- Output folder for plots ---
params.plot.folder = "../plots/part1";

if ~exist(params.plot.folder, "dir")
    mkdir(params.plot.folder);
end

%% --- Common setup ---
params.fileName = 'veldata.txt';

% Sampling frequency is 20 kHz
params.samplingFreq = 20000;

% Plot parameters
params.plot.fontSize = 14;
params.plot.legendFontSize = 12;
params.plot.insetFontSize = 10;
params.plot.lw = 1.5;
params.plot.lwsmall = 1.2;

% Load data
[u, U] = load_data(params.fileName, Inf);

% Mean velocity
params.U = U;

%% --- Frozen flow hypothesis ---
% Create the time vector based on the number of samples and the sampling frequency
time_vector = (0:length(u)-1)' / params.samplingFreq;

% Apply Taylor's frozen flow hypothesis: t = -x/U -> x = -t*U
% x values
params.x = -time_vector * params.U;

% Length 
params.L = abs(params.x(end));


%% --- Plot A ---
plotA = figure();
hold on; grid on;

% Plot total velocity u_tot(x) = u + U against x
plot(params.x, u + params.U, 'LineWidth', params.plot.lwsmall);

yline(params.U, '--r', 'Mean Velocity $U$', ...
     'Interpreter', 'latex', 'fontsize', params.plot.fontSize, ...
     'LineWidth', params.plot.lw, ...
     'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off', ...
     'LabelHorizontalAlignment', 'left');

xlabel('$x\;[m]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

ylabel('$u_\mathrm{tot}(x)\;[m/s]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

% Save the plot automatically
exportgraphics(plotA, fullfile(params.plot.folder, "plotA.png"), ...
    "Resolution", 300);

%% --- Turbulence intensity ---
params.I = sqrt(mean(u.^2)) / params.U;
fprintf('Mean velocity: %e\n', params.U);
fprintf('Turbulence intensity: %e\n', params.I);