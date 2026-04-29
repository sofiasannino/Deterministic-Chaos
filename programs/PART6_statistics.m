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
warning('Set the length increments and other relvant parameters');


%% --- Loop over dataset ---
plotE2 = figure();
plotE1 = figure();

params.pdf = cell(1, length(lincr));
params.pdf_gauss = cell(1, length(lincr));
params.pdf_incr = cell(1, length(lincr));
histNbins = 100;

for scale = 1:2  % 1: semilogy, 2: linear
    if (scale == 2)
        figure(plotE2);
    end
    for i = 1:length(lincr)
        s = subplot(2, 2, i);

        [x, d] = increment(u, lincr(i), params.samplingFreq, U);
        h = histogram(d, i * histNbins, 'Normalization', 'probability');
        increments = zeros(h.NumBins, 1);
        values = zeros(h.NumBins, 1);
        dx = h.BinWidth;

        for j = 1:h.NumBins
            increments(j) = (h.BinEdges(j)+h.BinEdges(j+1))/2;
            values(j) = h.Values(j);
        end

        values = values / dx;
        cla(s)

        if scale == 1
            % semilogy(..., 'LineWidth', params.plot.lw);
            ylim([1e-7 5e1])
        else
            % plot(..., 'LineWidth', params.plot.lw);
        end

        % bound = ...;
        % xlim([-bound, bound]);

        G = fit_gaussian(increments, d);
        hold on; grid on;
        % plot(..., 'LineWidth', params.plot.lw);
        % params.pdf_gauss{i} = G;   % Suggestion
        % params.pdf{i} = values;
        % params.pdf_incr{i} = increments;

        % Axis labels
        if i == 3 || i == 4
            xlabel('$\delta u_{||}\;[m/s]$', 'Interpreter', 'latex', ...
                'fontsize', params.plot.fontSize);
        end
        if i == 1 || i == 3
            ylabel('$P(\delta u_{||})$', 'Interpreter', 'latex', ...
                'fontsize', params.plot.fontSize);
        end

        % Legend
        if i == 4
            legend('Data PDF', 'Gaussian', 'Interpreter', 'latex', ...
                'fontsize', params.plot.legendFontSize); 
        end

        title(['$l=', num2str(lincr(i)), '\;m$'], 'Interpreter', ...
            'latex', 'fontsize', params.plot.fontSize);
    end

end
