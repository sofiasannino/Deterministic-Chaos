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

% l = ...
N = length(u);
L = l(N);
params.N = N;
params.L = L;


%% --- Energy spectrum ---
spectrumNorm = 1;
Ek = spectral_energy_density(u) / spectrumNorm;

% Define E(k), k > 0
% Ek = ... 

k = 1;

% Signal filtering
Ek_smooth = Ek;

warning('Make sure to fill in these variables and smooth the spectrum');


%% Plot C
plotC = figure();
% loglog(..., 'DisplayName', '$E(k)$ (smoothed)', ...
%     'LineWidth', params.plot.lw);
hold on; grid on;

% Add K41 scaling predictions
% ...

% Add estimates for integral and Kolmogorov length scales
params.eta_e = 1;
params.Lint_e = 1;

keta_e = 1;
kint_e = 1;
kLc = 1;

% Labels
xlabel('$k \;[m^{-1}]$', 'Interpreter', 'latex')
ylabel('$E(k) \; [m^3/s^2]$', 'Interpreter', 'latex')
legend('Interpreter', 'latex', 'FontSize', params.plot.legendFontSize);

% Plot correlation length scales Lci
% xline(..., ...
%     'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
%     'LabelVerticalAlignment', 'bottom', ...
%     'LabelHorizontalAlignment', 'left','HandleVisibility', 'off');



%% Parsival's theorem
warning('Make sure to check Parsival''s theorem');
