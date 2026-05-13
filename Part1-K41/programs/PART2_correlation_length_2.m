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
params.dx = U./params.samplingFreq; % dx translation in space for every sample: dx = tau * U = U/f
params.x = (0:length(u)-1)' * params.dx; % define x coordinates
params.L = params.x(end); % length of the interval on axis x

%% --- Correlation length ---
dl = params.dx;

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
lmax_plot = 1;
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
% compute integral length with trapezoid quadrature function
L_int = trapz(l(1:lmax_idx), Cl(1:lmax_idx)); 
fprintf('Integral scale: %e \n', L_int)
fprintf('Error: %e', norm(L_int - params.Lc)/L_int)