%% --- Move to this file's folder ---
cd(fileparts(mfilename('fullpath')));

%% --- Clear workspace ---
clc;
addpath("../functions/");

%% --- Output folder for plots ---
params.plot.folder = "../plots";

if ~exist(params.plot.folder, "dir")
    mkdir(params.plot.folder);
end
%% --- Common setup ---
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

% Here eta and lambda are not exact wavelengths of single Fourier modes,
% they are characteristic turbulent length scales obtained from physical
% estimates

% When we mark them on the spectrum, we only want to indicate the approximate
% wavenumber range where scales of size l become relevant. For this reason
% we use the standard turbulence convention k*l ~ 1, i.e. k = 1/l

% Using k = 2*pi/l would instead mark the Fourier mode whose wavelength is
% exactly l, which is not the interpretation used here
klambda = 2 * pi / params.lambda;
keta = 1 / params.eta;

fprintf('Dissipation rate (epsilon): %e\n', params.epsilon);
fprintf('Taylor length scale (lambda): %e\n', params.lambda);
fprintf('Taylor Reynolds number (Re_lambda): %.2f\n', params.Re_lambda);
fprintf('Outer Reynolds number (Re_outer): %.2f\n', params.Re_outer);
fprintf('Kolmogorov length scale (eta): %e\n', params.eta);
fprintf('Kolmogorov Reynolds number (Re_eta): %.2f\n', params.Re_eta);

%% --- Add lambda and eta scales to Plot C ---
plotC_fig = fullfile(params.plot.folder, "plotC.fig");
plotC_png = fullfile(params.plot.folder, "plotC_with_lambda_eta.png");

% Look for an already open Plot C
figs = findall(0, "Type", "figure", "Name", "PlotC");

if ~isempty(figs)
    plotC = figs(1);
else
    if exist(plotC_fig, "file")
        plotC = openfig(plotC_fig, "visible");
        plotC.Name = "PlotC";
    else
        warning("Plot C is neither open nor saved as plots/plotC.fig. Run Section 1.3 first.");
        plotC = [];
    end
end

if ~isempty(plotC)
    figure(plotC);

    % Get axes
    ax_all = findall(plotC, "Type", "axes");
    ax = ax_all(1);

    hold(ax, "on");

    % Avoid duplicated lines if the script is run multiple times
    delete(findall(ax, "Tag", "k_lambda_part4"));
    delete(findall(ax, "Tag", "k_eta_part4"));

    xline(ax, klambda, '--m', '$k_{\lambda}=1/\lambda$', ...
    'Interpreter', 'latex', ...
    'FontSize', params.plot.insetFontSize, ...
    'Color', 'black', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'right', ...
    'HandleVisibility', 'off', ...
    'Tag', 'k_lambda_part4');

    xline(ax, keta, '--c', '$k_{\eta}=1/\eta$', ...
    'Interpreter', 'latex', ...
    'FontSize', params.plot.insetFontSize, ...
    'Color', '#FF8C00', ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'right', ...
    'HandleVisibility', 'off', ...
    'Tag', 'k_eta_part4');

    % Save updated Plot C
    savefig(plotC, plotC_fig);
    exportgraphics(plotC, plotC_png, "Resolution", 300);
end