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

% l_S2  = ...
% l_S3  = ...

params.S2_exp  = 1;   % scaling exponents
params.S3_exp  = 1;
S2_exp_str = '1';
S3_exp_str = '1';

warning('Make sure to fill in the relevant values!');

%% --- Compute S2 and S3 ---
S2 = structure_function(u, l_S2, params.samplingFreq, U, 2);
S3 = structure_function(u, l_S3, params.samplingFreq, U, 3);


%% --- Smooth data ---
S2_smooth = S2;
S3_smooth = S3;
warning('Make sure to use appropriate smoothing');


%% --- Plot S2 ---
plotF = figure();

% loglog(..., 'LineWidth', params.plot.lw);
hold on; grid on;
% loglog(..., 'LineWidth', ...
%     params.plot.lwsmall);

xlabel('$l\;[m]$', 'Interpreter', 'latex', ...
    'fontsize', params.plot.fontSize);
ylabel('$S_2(l)\;[m^2/s^2]$', 'Interpreter', 'latex', ...
    'fontsize', params.plot.fontSize);
% legend(..., ...
%     'Interpreter', 'latex', ...
%     'fontsize', params.plot.legendFontSize, ...
%     'Location','NorthWest');

% xline(...,
%     'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
%     'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');



%% --- Plot S3 ---
plotG = figure();
S3_magnitude = 1;

% loglog( ..., 'LineWidth', params.plot.lw);
hold on; grid on;

% loglog(..., 'LineWidth', ...
    % params.plot.lwsmall);

xlabel('$l\;[m]$', 'Interpreter', 'latex', ...
    'fontsize', params.plot.fontSize);
ylabel('$S_3(l)\;[m^3/s^3]$', 'Interpreter', 'latex', ...
    'fontsize', params.plot.fontSize);
% legend(, ...
%     'Interpreter', 'latex', ...
%     'fontsize', params.plot.legendFontSize, ...
%     'Location','NorthWest');

% xline(...
%     'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
%     'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');



%% --- Estimate ε from S2 and S3 ---
% Estimate from S3
epsilon_S3 = 1;
    
% Estimate from S2
C2 = 1;
epsilon_S2 = 1;

warning('Make sure to use fill those parameters');


%% --- Add result on plots F and G ---
figure(plotF);
% loglog(...,    'LineWidth', params.plot.lw, 'DisplayName', ...
%     '$y = C_2 \left( \langle \epsilon \rangle l \right)^{2/3}$');

figure(plotG);
% loglog(...
%     'LineWidth', params.plot.lw, 'DisplayName', ...
%     '$y = \frac{4}{5} \langle \epsilon \rangle l $');

warning('Make sure to use fill those plotting lines');