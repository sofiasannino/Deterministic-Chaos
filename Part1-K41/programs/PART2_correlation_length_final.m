%% --- Clear workspace ---
clc;
clear;
addpath("../functions/");


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
params.L = params.x(end); % length of the interval on axis x

%% --- Correlation length ---
dl = abs(params.dx); %absolute displacement

lmax_autocorr = 1000;
Cl = autocorrelation(u, dl, lmax_autocorr);


% Length vectors
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
err_zero1 = norm(L_zero1 - params.Lc) / L_zero1;
err_min   = norm(L_min   - params.Lc) / L_min;
err_zero2 = norm(L_zero2 - params.Lc) / L_zero2;

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