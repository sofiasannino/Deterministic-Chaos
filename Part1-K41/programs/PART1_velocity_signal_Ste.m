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

% Apply Taylor's frozen flow hypothesis: t = -x/U -> x = -t*U
params.x = -time_vector * params.U;
params.L = abs(params.x(end));

%% --- Plot A ---
plotA = figure();
hold on; grid on;

% Plot total velocity u_tot(x) against x
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


%% --- Turbulence intensity ---
params.I = sqrt(mean(u.^2)) / abs(params.U);
fprintf('Mean velocity U: %e\n', params.U);
fprintf('Mean of velocity fluctuations u: %e\n', mean(u));
fprintf('Turbulence intensity I: %e\n', params.I);

%% ADDED --- Median deviation of local averages from zero ---
u0 = u(:);                  % make column vector
N = length(u);
f=20000;

% --- values of n ---

% Smallest window: 1 sample
n_min = 1;

% Largest window: avoid using windows comparable to the full signal length
n_max = floor(N/2);

% Use logarithmically spaced window sizes
% n_values = unique(round(logspace(log10(n_min), log10(n_max), 100)));
n_values=unique(floor(linspace(10,n_max,150)));
% Preallocate
mean_abs_median = zeros(size(n_values));
time_windows = n_values / f;              % [s]
ell_values = abs(U) * time_windows;       % [m], using Taylor hypothesis

% --- Compute madian of absolute moving average for each n ---

for i = 1:length(n_values)

    n = n_values(i);

    % Moving average over the current point and the previous n-1 points
    local_mean = movmean(u0, [n-1 0]);

    % Remove the first n-1 values, where the window is not complete
    local_mean = local_mean(n:end);

    % Median deviation from zero
    mean_abs_median(i) = median(abs(local_mean));

end

tol = 1e-2;

acceptable_idx = find(mean_abs_median < tol, 1, 'first');
acceptable_ell = ell_values(acceptable_idx);
acceptable_n = n_values(acceptable_idx);
acceptable_error = mean_abs_median(acceptable_idx);
% --- Plot versus spatial length ---
figure;
ax1 = axes;
semilogy(ax1, ell_values, mean_abs_median, 'LineWidth', 1.5);
hold(ax1, 'on');
grid(ax1, 'on');
box(ax1, 'on');

xlabel(ax1, '$\ell = |U| n/f$ [m]', 'Interpreter', 'latex');
ylabel(ax1, '$\mathrm{median}\left(|\overline{u}_n|\right)$ [m/s]', ...
       'Interpreter', 'latex');

set(ax1, 'FontSize', 14);

% Add threshold and first acceptable length
yline(ax1, tol, ':', '$10^{-2}$', ...
      'Interpreter', 'latex', ...
      'LabelHorizontalAlignment', 'left');

if ~isempty(acceptable_idx)

    xline(ax1, acceptable_ell, '--', ...
          sprintf('$\\ell_* = %.3g$ m', acceptable_ell), ...
          'Interpreter', 'latex', ...
          'LabelOrientation', 'horizontal', ...
          'LabelVerticalAlignment', 'bottom');

    plot(ax1, acceptable_ell, acceptable_error, 'o', ...
         'MarkerSize', 7, ...
         'LineWidth', 1.5);

    text(ax1, acceptable_ell, acceptable_error, ...
         sprintf('  $n_* = %d$', acceptable_n), ...
         'Interpreter', 'latex', ...
         'VerticalAlignment', 'bottom');

end

% --- Add second x-axis showing number of samples n ---

% Leave more vertical space above the plot
ax1.Position = [0.13 0.14 0.775 0.68];

ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', ...
           'Color', 'none', ...
           'YTick', [], ...
           'XScale', ax1.XScale, ...
           'XLim', ax1.XLim);

% Choose ticks on the bottom axis
top_ticks_l = ax1.XTick;

% Convert spatial length back to number of samples:
% ell = |U| n/f  =>  n = ell f / |U|
top_ticks_n = top_ticks_l * f / abs(U);

ax2.XTick = top_ticks_l;
ax2.XTickLabel = compose('%.0f', top_ticks_n);

xlabel(ax2, 'Window size $n$ [samples]', 'Interpreter', 'latex');

% Move the top-axis label slightly upward
ax2.XLabel.Units = 'normalized';
ax2.XLabel.Position(2) = 1.08;

% Main title, moved above the second x-axis label
t = title(ax1, 'Median deviation from zero of local averages', ...
          'Interpreter', 'latex');

t.Units = 'normalized';
t.Position(2) = 1.17;

% Make sure axes stay aligned
linkaxes([ax1 ax2], 'x');
