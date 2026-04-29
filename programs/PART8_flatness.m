%% --- Clear workspace ---
clc;
clear;
addpath("../functions/");


%% --- Common setup ---
% Define any quantity you need for this part in here
% l_S2 = ...
% l_S3 = ...


%% --- Compute flatness ---
S4 = structure_function(u, l_S2, params.samplingFreq, U, 4);
S6 = structure_function(u, l_S2, params.samplingFreq, U, 6);


%% --- Flatness and Hyper-flatness ---
S4_smooth = S4;
S6_smooth = S6;
% F = ...;
% H = ...;


%% --- Plot Flatness ---
plotH = figure();
% loglog(..., 'LineWidth', params.plot.lw);
hold on; grid on;

% loglog(..., ...
    % 'LineWidth', params.plot.lwsmall);

% yline( ...
    % 'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
    % 'LabelHorizontalAlignment','left');


xlabel('$l\;[m]$', 'Interpreter', 'latex', 'fontsize', ...
    params.plot.fontSize);
ylabel('$F(l) = \frac{S_4(l)}{S_2(l)^2}$', 'Interpreter', ...
    'latex', 'fontsize', params.plot.fontSize);
% legend(..., ...
%     'Interpreter', 'latex', 'fontsize', params.plot.legendFontSize);

% xline(..., ...
    % 'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
    % 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');


%% --- Correction to energy spectrum ---
figure(plotC);
cmap = colororder();  % Get the default colors if needed

% plot(..., ...
%     'LineWidth', params.plot.lwsmall, ...
%     'DisplayName', '$\propto k^{-51/30}$ ($\beta$-model)', ...
%     'Color', cmap(2, :));
% plot(..., ...
%     'LineWidth', params.plot.lwsmall, ...
%     'DisplayName', '$\propto k^{-53/30}$ (bifrac.)', ...
%     'Color', cmap(3, :));


%% --- Plot Hyper-flatness ---
plotI = figure();
% loglog(..., 'LineWidth', params.plot.lw);
hold on; grid on;

% loglog(..., ...
%     'LineWidth', params.plot.lwsmall);
% loglog(..., ...
%     'LineWidth', params.plot.lwsmall);

% yline(..., ...
%     'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
%     'LabelHorizontalAlignment','left');

xlabel('$l\;[m]$', 'Interpreter', 'latex', 'fontsize', ...
    params.plot.fontSize);
ylabel('$H(l) = \frac{S_6(l)}{S_2(l)^3}$', 'Interpreter', 'latex', ...
    'fontsize', params.plot.fontSize);
% legend(..., ...
%     'Interpreter', 'latex', 'fontsize', params.plot.legendFontSize);

% xline..., ...
%     'Interpreter', 'latex', 'fontsize', params.plot.insetFontSize, ...
%     'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

