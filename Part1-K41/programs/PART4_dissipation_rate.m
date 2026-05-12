%% --- Clear workspace ---
clc;
addpath("../functions/");


%% --- Common setup ---
% Define any quantity you need for this part in here
params.fileName = 'veldata.txt';
params.Lc = 0.3573;
params.plot.fontSize = 14;
params.plot.legendFontSize = 12;
params.plot.insetFontSize = 10;
params.plot.lw = 1.5;
params.plot.lwsmall = 1.2;

[u, U] = load_data(params.fileName, Inf);


%% --- Dissipation and Reynolds numbers ---
% Fluid's viscosity
nu = 1.51e-5;
L_outer = params.Lc;

u_rms = sqrt(mean(u.^2));

params.epsilon = 0.5 * (u_rms^3) / L_outer;
params.lambda = sqrt(15 * nu * (u_rms^2) / params.epsilon);
params.eta = (nu^3 / params.epsilon)^(1/4);
params.Re_lambda = (u_rms * params.lambda) / nu;
params.Re_outer = (u_rms * L_outer) / nu;

% By definition this is always exactly 1, but we calculate it   
u_eta = (params.epsilon * nu)^(1/4);
params.Re_eta = (u_eta * params.eta) / nu;

klambda = 2 * pi / params.lambda;
keta = 2 * pi / params.eta;

fprintf('Dissipation rate (epsilon): %e m^2/s^3\n', params.epsilon);
fprintf('Taylor length scale (lambda): %e m\n', params.lambda);
fprintf('Taylor Reynolds number (Re_lambda): %.2f\n', params.Re_lambda);
fprintf('Outer Reynolds number (Re_outer): %.2f\n', params.Re_outer);
fprintf('Kolmogorov Reynolds number (Re_eta): %.2f\n', params.Re_eta);

% Add scales to plot C
figure(plotC)
try
    figure(plotC)
    hold on;
    
    xline(klambda, '--b', '$k_{\lambda}$', ...
        'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
        'LabelVerticalAlignment', 'bottom', ...
        'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');
catch
    warning('plotC is not open. Run Section 1.3 first to add these lines to the spectrum.');
end