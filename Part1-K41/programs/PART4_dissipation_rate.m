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
params.plot.fontSize = 20;
params.plot.legendFontSize = 18;
params.plot.insetFontSize = 16;
params.plot.lw = 2.1;
params.plot.lwsmall = 1.8;

[u, U] = load_data(params.fileName, Inf);

%% --- Dissipation and Reynolds numbers ---
% Fluid's viscosity
nu = 1.51e-5;
L_outer = 0.3573;
u_rms = sqrt(mean(u.^2));

params.epsilon = 0.5 * (u_rms^3) / L_outer;
params.lambda = sqrt(15 * nu * (u_rms^2) / params.epsilon);
params.eta = (nu^3 / params.epsilon)^(1/4);
params.Re_lambda = (u_rms * params.lambda) / nu;
params.Re_outer = (u_rms * L_outer) / nu;

% By definition this is always exactly 1, but we calculate it   
u_eta = (params.epsilon * nu)^(1/4);
params.Re_eta = (u_eta * params.eta) / nu;

%% --- Add lambda and eta scales to Plot C ---
klambda = 2 * pi / params.lambda;

keta = 1 / params.eta;
keta_min = pi / params.eta / 4;
keta_max = pi / params.eta / 30;

% Look for an already open plot 1C
figs = findall(0, "Type", "figure", "Name", "Plot1C");

if ~isempty(figs)
    plotC = figs(1);
else
    if exist(fullfile(params.plot.folder, "plot1C.fig"), "file")
        plot1C = openfig(fullfile(params.plot.folder, "plot1C.fig"), "visible");
        plot1C.Name = "Plot1C";
    else
        warning("Plot 1C is neither open nor saved as plots/plot1C.fig. Run Section 1.3 first.");
        plotC = [];
    end
end

if ~isempty(plot1C)
    figure(plot1C);

    % Get axes
    ax_all = findall(plot1C, "Type", "axes");
    ax = ax_all(1);
    hold(ax, "on");

    % Avoid duplicated lines if the script is run multiple times
    delete(findall(ax, "Tag", "k_lambda_part4"));
    delete(findall(ax, "Tag", "k_eta_part4"));
    delete(findall(ax, "Tag", "k_pope_part4"));
    delete(findall(ax, "Tag", "k_text_part4"));

    xline(ax, klambda, '--m', '$k_{\lambda}$', ...
    'Interpreter', 'latex', ...
    'FontSize', params.plot.insetFontSize, ...
    'Color', 'black', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off', ...
    'Tag',"k_lambda_part4" );

    xline(ax, keta, '--', '$k_{\eta}$', ...
    'Interpreter', 'latex', ...
    'FontSize', params.plot.insetFontSize, ...
    'Color', '#FF8C00', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off', ...
    'Tag', 'k_eta_part4');

    % Find the existing legend in the current axes
    lgd = findobj(plot1C, 'Type', 'Legend');
    if ~isempty(lgd)
        lgd.Location = 'southwest'; 
    end

    p_eta_range = xregion(keta_min, keta_max, ...
    'FaceColor', [1.00 0.78 0.25], ...
    'FaceAlpha', 0.45, ...
    'EdgeColor', 'none', ...
    'HandleVisibility', 'off', ...
    'Tag', 'k_pope_part4');

    uistack(p_eta_range, 'bottom');

    yl_F = ylim;
    text(sqrt(keta_min*keta_max), yl_F(1)*(yl_F(2)/yl_F(1))^0.08, ...
        '$\eta\cdot8 < l < \eta\cdot60$', ...
        'Interpreter', 'latex', ...
        'fontsize', params.plot.insetFontSize, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'Tag', 'k_text_part4');

    % Save the plot automatically
    if isgraphics(plot1C)
        savefig(plot1C, fullfile(params.plot.folder, "plot1C.fig"));
        exportgraphics(plot1C, fullfile(params.plot.folder, "plot1C.png"), "Resolution", 300);
    else
        disp('Skipped exporting plot1C (figure was closed).');
    end
end

fprintf('================ Part 1.4 results ================\n');
fprintf('Dissipation rate (epsilon): %e\n', params.epsilon);
fprintf('Taylor length scale (lambda): %e\n', params.lambda);
fprintf('Taylor Reynolds number (Re_lambda): %.2f\n', params.Re_lambda);
fprintf('Outer Reynolds number (Re_outer): %.2f\n', params.Re_outer);
fprintf('Kolmogorov length scale (eta): %e\n', params.eta);
fprintf('Kolmogorov Reynolds number (Re_eta): %.2f\n', params.Re_eta);