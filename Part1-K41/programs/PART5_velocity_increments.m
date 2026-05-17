%% --- Move to this file's folder ---
cd(fileparts(mfilename('fullpath')));

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
params.plot.fontSize = 14;
params.plot.legendFontSize = 12;
params.plot.insetFontSize = 10;
params.plot.lw = 1.5;
params.plot.lwsmall = 1.2;

params.fileName = 'veldata.txt'; %
params.samplingFreq = 20e3;  %

[u, U] = load_data(params.fileName, Inf);

lincr = [1e-3, 1e-2, 1e-1, 10]; %
lincr_labels = {'1 mm', '1 cm', '10 cm', '10 m'};

%% --- Loop through dataset ---
plotD = figure();
maxx_plot = 1600;

for i = 1:length(lincr)  % loop over l values     
    [x, d] = increment(u, lincr(i), params.samplingFreq, U);

    % With the convention x = -Ut, we display the interval [-maxx_plot, 0].
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
        xlabel('$x\;[m]$', 'Interpreter', 'latex', 'fontsize', ...
            params.plot.fontSize);
    end

    if i == 1 || i == 3
        ylabel('$\delta u_{||}\;[m/s]$', 'Interpreter', 'latex', ...
            'fontsize', params.plot.fontSize);
    end

    % Title
    subplot_title = ['$l=', num2str(lincr(i)), '\;m$'];
    title(subplot_title, 'Interpreter', 'latex', 'fontsize', ...
        params.plot.legendFontSize);

    % print mean increment (should be ~0)
    disp(['l = ', num2str(lincr(i)), ': mean = ', num2str(mean(d))]);
end

exportgraphics(plotD, fullfile(params.plot.folder, "plotD.png"), ...
    "Resolution", 300);