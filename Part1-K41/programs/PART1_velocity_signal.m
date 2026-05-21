%% --- Clear workspace ---
clc;
clear;

%% --- Robust Path Setup ---
% Find the directory of this script, whether running via F5 or F9
if isempty(mfilename)
    % Fallback for when running code line-by-line (F9) in the Editor
    scriptDir = fileparts(matlab.desktop.editor.getActiveFilename);
else
    % Standard execution (F5 or command line)
    scriptDir = fileparts(mfilename('fullpath'));
end

% Build absolute paths relative to the script's directory
functionsDir = fullfile(scriptDir, '..', 'functions');
dataDir = fullfile(scriptDir, '..', 'data');
params.plot.folder = fullfile(scriptDir, '..', 'plots');

% Add dependencies to path
addpath(functionsDir);
addpath(dataDir);

% Create the plot folder if it doesn't exist
if ~exist(params.plot.folder, 'dir')
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

% Load data
[u, U] = load_data(params.fileName, Inf);

% Mean velocity
params.U = U;

%% --- Frozen flow hypothesis ---
% Create the time vector based on the number of samples and the sampling frequency
time_vector = (0:length(u)-1)' / params.samplingFreq;

% Apply Taylor's frozen flow hypothesis: t = -x/U -> x = -t*U
params.x = -time_vector * params.U;
params.L = abs(params.x(end));

%% --- Plot 1A ---
plot1A = figure("Name", "Plot1A");
hold on; grid on;

% Plot total velocity u_tot(x) = u + U against x
plot(params.x, u + params.U, 'LineWidth', params.plot.lwsmall);

yline(params.U, '--r', 'Mean Velocity $U$', ...
     'Interpreter', 'latex', 'fontsize', params.plot.fontSize, ...
     'LineWidth', params.plot.lw, ...
     'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off', ...
     'LabelHorizontalAlignment', 'left');

xlabel('$x\;[\mathrm{m}]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

ylabel('$u_\mathrm{tot}(x)\;[\mathrm{m/s}]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);

% Save the plot automatically
if isgraphics(plot1A)
    exportgraphics(plot1A, fullfile(params.plot.folder, "plot1A.png"), "Resolution", 300);
else
    disp('Skipped exporting plot1A (figure was closed).');
end

%% --- Turbulence intensity ---
params.I = sqrt(mean(u.^2)) / params.U;

%% --- Display results ---
fprintf('================ Part 1.1 results ================\n');
fprintf('Mean velocity: %e\n', params.U);
fprintf('Turbulence intensity: %e\n', params.I);