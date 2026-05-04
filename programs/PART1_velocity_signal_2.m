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


%% --- Frozen flow hypothesis ---
params.U = U;
params.x = U./params.samplingFreq; % dx translation in space for every sample: dx = tau * U
% define x coordinates
xx = (0:length(u)-1)' * params.x;
params.L = xx(end); % length of the interval on axis x

% define total velocity
u_tot = u + U.*ones(length(u), 1);

%% --- Plot A ---
plotA = figure();
hold on; grid on;
plot(xx, u_tot, 'LineWidth', 1);
%xlim([xx(1), xx(end)]); % to center the figure 

yline(U, '--r', 'Mean Velocity $U$',...
    'Interpreter', 'latex', 'fontsize', params.plot.fontSize, ...
     'LineWidth', params.plot.lw, ...
     'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off', ...
     'LabelHorizontalAlignment', 'left');

xlabel('$x\;[m]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);
ylabel('$u_\mathrm{tot}(x)\;[m/s]$', 'Interpreter', 'latex', ...
    'FontSize', params.plot.fontSize);


%% --- Turbulence intensity ---
u_squared_avg = mean(u.^2);
params.I = sqrt(u_squared_avg) ./ U;
fprintf('Turbulence intensity: %e\n', params.I);