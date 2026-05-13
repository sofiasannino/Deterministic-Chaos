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

params.fileName = 'veldata.txt'; %
params.samplingFreq = 20e3; %
[u, U] = load_data(params.fileName, Inf);

lincr = [1e-3, 1e-2, 1e-1, 10]; %


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
        h_exp = histogram(d, i * histNbins, 'Normalization', 'probability');
        increments = zeros(h_exp.NumBins, 1);
        values = zeros(h_exp.NumBins, 1);
        dx = h_exp.BinWidth;

        for j = 1:h_exp.NumBins
            increments(j) = (h_exp.BinEdges(j)+h_exp.BinEdges(j+1))/2;
            values(j) = h_exp.Values(j);
        end

        values = values / dx;
        cla(s)

        if scale == 1
            semilogy(increments, values, 'LineWidth', params.plot.lw);
            ylim([1e-7 5e1])
        else
            plot(increments, values, 'LineWidth', params.plot.lw);
        end

        bound = 6*std(d); %
        xlim([-bound, bound]);

        G = fit_gaussian(increments, d);
        hold on; grid on;
        plot(increments, G, 'LineWidth', params.plot.lw);
        params.pdf_gauss{i} = G;   % Suggestion
        params.pdf{i} = values;
        params.pdf_incr{i} = increments;

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

%% Inspecting if H2 stands 
h_exp = zeros(length(lincr)-1, 5);

for i = 1:(length(lincr)-1)

    du_i   = params.pdf_incr{i};
    du_i_1 = params.pdf_incr{i+1};

    pdf_i   = params.pdf{i};
    pdf_i_1 = params.pdf{i+1};

    % common support for interpolation
    common = du_i_1 >= min(du_i) & du_i_1 <= max(du_i);

    du_common = du_i_1(common);
    pdf2 = pdf_i_1(common);

    % interpolate pdf_i in the grid of pdf_i_1
    pdf1 = interp1(du_i, pdf_i, du_common, 'spline');

    % remove NaNs etc
    valid = isfinite(pdf1) & isfinite(pdf2) & pdf1 > 0 & pdf2 > 0;

    du_valid = du_common(valid);
    pdf1_valid = pdf1(valid);
    pdf2_valid = pdf2(valid);

    for p = 1:5

        % moments
        m1 = trapz(du_valid, abs(du_valid).^p .* pdf1_valid);
        m2 = trapz(du_valid, abs(du_valid).^p .* pdf2_valid);

        % estimate h 
        h_exp(i,p) = (1/p) * log10(m2/m1) / log10(lincr(i+1)/lincr(i));

    end
end
