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


%% --- Dissipation and Reynolds numbers ---
nu = 1;                         % Fluid's viscosity
L_outer = 1;

params.epsilon = 1;
params.lambda = 1;
params.Re_lambda = 1;
params.Re_outer = 1;

klambda = 1;
keta = 1;


% Add scales to plot C
figure(plotC)
% xline(..., ...
%     'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
%     'LabelVerticalAlignment', 'bottom', ...
%     'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');
%     % you might want to adjust the horizontal label alignement


warning('Make sure that you filled all of the above numbers')