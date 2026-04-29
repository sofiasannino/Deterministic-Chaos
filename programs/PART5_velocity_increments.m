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

% lincr = [...];
params.samplingFreq = 1;

warning('Set the length increments and other relvant parameters');


%% --- Loop through dataset ---
plotD = figure();
maxx_plot = 1600;

for i = 1:length(lincr)  % loop over l values     
    [x, d] = increment(u, lincr(i), params.samplingFreq, U);
    x = x(x <= maxx_plot);
    subplot(2, 2, i)
    % plot(...);
    miny = min(d);
    maxy = max(d);
    if abs(miny) > maxy
        maxy = abs(miny);
    elseif maxy > abs(miny)
        miny = -maxy;
    end
    % xlim([...]);
    % ylim([...]);
    % yline(0, '-', ...);

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
