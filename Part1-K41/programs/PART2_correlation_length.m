%% --- Clear workspace ---
clc;
clear;
addpath("../functions/");

%% --- Output folder for plots ---
params.plot.folder = "plots";

if ~exist(params.plot.folder, "dir")
    mkdir(params.plot.folder);
end
%% --- Common setup ---
% Define any quantity you need for this part in here
params.fileName = 'veldata.txt';
params.samplingFreq = 20000; 
params.plot.fontSize = 14;
params.plot.legendFontSize = 12;
params.plot.insetFontSize = 10;
params.plot.lw = 1.5;
params.plot.lwsmall = 1.2;

[u, U] = load_data(params.fileName, Inf);
params.U = U;
params.dx = - U./params.samplingFreq; % signed Taylor displacement
params.x = (0:length(u)-1)' * params.dx; % define x coordinates
params.L = abs(params.x(end)); %absolute length of the interval on axis x

%% --- Correlation length ---
dl = abs(params.dx); %absolute displacement

lmax_autocorr = 1000;
Cl = autocorrelation(u, dl, lmax_autocorr); %this function calculates correlation in accordance with our convension: a correlation of +l means we take all couples of velocities with index distance freq * l / U and they are computed excluding the 'physical right tail' of the velocities, meaning that for an index q, velocities 0,1,..., q-1 (which physically are on the right of our domain, close to zero) won't have be considered, correctly so since we don't have data for velocities at physical distance +l from them. 


% Absolute length vectors
l = 0:dl:(dl*(length(Cl) - 1));

% Find intersections
critic_value = 1/exp(1);
Lc_idx = find(Cl <= critic_value, 1, 'first');
params.Lc = l(Lc_idx);
fprintf('critical length: %e \n', params.Lc)


%% --- Plot B ---
lmax_plot = 10;
lmax_idx = find(l >= lmax_plot, 1);
if isempty(lmax_idx)
    lmax_idx = length(l);
end

plotB = figure();
hold on; grid on;
plot(l(1:lmax_idx), Cl(1:lmax_idx), 'LineWidth', params.plot.lw);

xlabel('$l\;[m]$', 'Interpreter', 'latex', 'FontSize', params.plot.fontSize);
ylabel('$C(l)$', 'Interpreter', 'latex', 'FontSize', params.plot.fontSize);

% Add Lc to plots
xline(params.Lc, '--r', 'Correlation length Lc', ...
    'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% Save plot automatically
exportgraphics(plotB, fullfile(params.plot.folder, "plotB.png"), ...
    "Resolution", 300);
%% --- Integral scale ---
% First zero crossing: Cl becomes negative
idx_zero1 = find(Cl <= 0, 1, 'first');

% Second zero crossing: Cl becomes positive again
idx_zero2 = idx_zero1 + find(Cl(idx_zero1+1:end) >= 0, 1, 'first');


% Minimum of Cl between first and second zero crossing
[~, idx_min_rel] = min(Cl(idx_zero1:idx_zero2));
idx_min = idx_zero1 + idx_min_rel - 1;
fprintf('minimum value of Cl: %e \n', Cl(idx_min))

% Integrals
L_zero1 = trapz(l(1:idx_zero1), Cl(1:idx_zero1));
L_min   = trapz(l(1:idx_min),   Cl(1:idx_min));
L_zero2 = trapz(l(1:idx_zero2), Cl(1:idx_zero2));

% Errors
err_zero1 = abs(L_zero1 - params.Lc) / params.Lc;
err_min   = abs(L_min   - params.Lc) / params.Lc;
err_zero2 = abs(L_zero2 - params.Lc) / params.Lc;

% Print results
fprintf('\n================ Integral length scales ================\n')

fprintf('\n[Stop at first zero crossing]\n')
fprintf('L_int: %e\n', L_zero1)
fprintf('Error: %e\n', err_zero1)

fprintf('\n[Stop at minimum]\n')
fprintf('L_int: %e\n', L_min)
fprintf('Error: %e\n', err_min)

fprintf('\n[Stop at second zero crossing]\n')
fprintf('L_int: %e\n', L_zero2)
fprintf('Error: %e\n', err_zero2)

fprintf('\n========================================================\n')

fprintf('\nValue used in the results table:\n')
fprintf('Lc    = %e m\n', params.Lc)
fprintf('Lint  = %e m\n', L_zero1)
fprintf('Relative difference |Lint - Lc|/Lc = %e\n', err_zero1)