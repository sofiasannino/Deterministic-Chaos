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

[u, U] = load_data(params.fileName, Inf);

lincr = [1e-3, 1e-2, 1e-1, 10];
lincr_labels = {'1 mm', '1 cm', '10 cm', '10 m'};

%% --- Loop through dataset ---
plot1D = figure();
maxx_plot = 100;

mean_increments = zeros(length(lincr), 1);

for i = 1:length(lincr)   
    [x, d] = increment(u, lincr(i), params.samplingFreq, U);

    % With the convention x = -Ut, we display the interval [-maxx_plot, 0]
    plot_idx = (x >= -maxx_plot) & (x <= 0);
    x_plot = x(plot_idx);
    d_plot = d(plot_idx);

    subplot(2, 2, i)
    plot(x_plot, d_plot, 'LineWidth', params.plot.lwsmall);

    miny = min(d_plot);
    maxy = max(d_plot);

    if abs(miny) > maxy
        maxy = abs(miny);
    elseif maxy > abs(miny)
        miny = -maxy;
    end

    xlim([-maxx_plot 0]);
    ylim([miny maxy]);
    yline(0, '-', 'HandleVisibility', 'off');

    % Axes labels
    if i == 3 || i == 4
        xlabel('$x\;[\mathrm{m}]$', 'Interpreter', 'latex', 'fontsize', ...
            params.plot.fontSize);
    end

    if i == 1 || i == 3
        ylabel('$\delta u_{||}\;[\mathrm{m/s}]$', 'Interpreter', 'latex', ...
            'fontsize', params.plot.fontSize);
    end

    % Title
    subplot_title = ['$l=', num2str(lincr(i)), '\; \mathrm{m}$'];
    title(subplot_title, 'Interpreter', 'latex', 'fontsize', ...
        params.plot.legendFontSize);

    mean_increments(i) = mean(d);
end

% Save the plot automatically
if isgraphics(plot1D)
    exportgraphics(plot1D, fullfile(params.plot.folder, "plot1D.png"), "Resolution", 300);
else
    disp('Skipped exporting plot1D (figure was closed).');
end

%% --- Display results ---
fprintf('================ Part 1.5 results ================\n');

for i = 1:length(lincr)
    fprintf('l = %e, mean = %e\n', lincr(i), mean_increments(i));
end