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

[u, U] = load_data(params.fileName, Inf);
lincr = [1e-3, 1e-2, 1e-1, 10];

%% --- Loop over dataset ---
plot1Ev1 = figure("Name", "Plot1Ev1");
plot1Ev2 = figure("Name", "Plot1Ev2");

params.pdf = cell(1, length(lincr));
params.pdf_gauss = cell(1, length(lincr));
params.pdf_incr = cell(1, length(lincr));
histNbins = 100;

for scale = 1:2
    if scale == 1
        figure(plot1Ev1);
    else
        figure(plot1Ev2);
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
            semilogy(increments, values, 'LineWidth', params.plot.lw);
            ylim([1e-7 5e1])
        else
            plot(increments, values, 'LineWidth', params.plot.lw);
        end

        bound = 6*std(d);
        xlim([-bound, bound]);

        G = fit_gaussian(increments, d);
        hold on; grid on;
        plot(increments, G, 'LineWidth', params.plot.lw);
        params.pdf_gauss{i} = G;
        params.pdf{i} = values;
        params.pdf_incr{i} = increments;

        % Axis labels
        if i == 3 || i == 4
            xlabel('$\delta u_{||}\;[\mathrm{m/s}]$', 'Interpreter', 'latex', ...
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

        title(['$l=', num2str(lincr(i)), '\; \mathrm{m}$'], 'Interpreter', ...
            'latex', 'fontsize', params.plot.fontSize);
    end

end

% Save the plot automatically
if isgraphics(plot1Ev1)
    exportgraphics(plot1Ev1, fullfile(params.plot.folder, "plot1Ev1.png"), "Resolution", 300);
else
    disp('Skipped exporting plot1Ev1 (figure was closed).');
end

if isgraphics(plot1Ev2)
    exportgraphics(plot1Ev2, fullfile(params.plot.folder, "plot1Ev2.png"), "Resolution", 300);
else
    disp('Skipped exporting plot1Ev2 (figure was closed).');
end

%% Inspecting if H2 stands
h_exp = zeros(length(lincr)-1, 5);

for i = 1:(length(lincr)-1)
    du_i   = params.pdf_incr{i};
    du_i_1 = params.pdf_incr{i+1};

    pdf_i   = params.pdf{i};
    pdf_i_1 = params.pdf{i+1};

    % Common support for interpolation
    common = du_i_1 >= min(du_i) & du_i_1 <= max(du_i);

    du_common = du_i_1(common);
    pdf2 = pdf_i_1(common);

    % Interpolate pdf_i in the grid of pdf_i_1
    pdf1 = interp1(du_i, pdf_i, du_common, 'spline');

    % Remove NaNs etc
    valid = isfinite(pdf1) & isfinite(pdf2) & pdf1 > 0 & pdf2 > 0;

    du_valid = du_common(valid);
    pdf1_valid = pdf1(valid);
    pdf2_valid = pdf2(valid);

    for p = 1:5

        % Compute moments
        m1 = trapz(du_valid, abs(du_valid).^p .* pdf1_valid);
        m2 = trapz(du_valid, abs(du_valid).^p .* pdf2_valid);

        % Estimate h 
        h_exp(i,p) = (1/p) * log10(m2/m1) / log10(lincr(i+1)/lincr(i));
    end
end

plot1Extra2 = figure();
hold on; grid on;

p_vals = 1:5;
markers = {'-o', '-s', '-^'};

labels = {
    sprintf('$l$: %g m $\\rightarrow$ %g m', lincr(1), lincr(2)), ...
    sprintf('$l$: %g m $\\rightarrow$ %g m', lincr(2), lincr(3)), ...
    sprintf('$l$: %g m $\\rightarrow$ %g m', lincr(3), lincr(4))
};

% Plot the estimated h for each interval
for i = 1:size(h_exp, 1)
    plot(p_vals, h_exp(i,:), markers{i}, 'LineWidth', params.plot.lw, 'DisplayName', labels{i});
end

% Add the K41 theoretical prediction (h = 1/3)
yline(1/3, 'k--', 'K41 ($h = 1/3$)', 'LineWidth', params.plot.lw, ...
    'Interpreter', 'latex', 'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xlabel('Moment order $p$', 'Interpreter', 'latex', 'FontSize', params.plot.fontSize);
ylabel('Estimated scaling exponent $h$', 'Interpreter', 'latex', 'FontSize', params.plot.fontSize);
title('Testing K41 H2: Self-Similarity', 'Interpreter', 'latex', 'FontSize', params.plot.fontSize);
legend('Interpreter', 'latex', 'FontSize', params.plot.legendFontSize, 'Location', 'northeast');

% Save the plot automatically
if isgraphics(plot1Extra2)
    exportgraphics(plot1Extra2, fullfile(params.plot.folder, "plot1Extra2.png"), "Resolution", 300);
else
    disp('Skipped exporting plot1Extra2 (figure was closed).');
end