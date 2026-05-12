close all; clear; clc

%% inputs
addpath('../functions/')        % folder containing functions
nu = 6.5e-4;                    % hyperviscosity parameter
N = 128;                        % spatial resolution
t_study = 2500.0 * nu;          % analysis time period
dt = 0.1 * nu;                  % time step size for time integration
dt_store = 1.0 * nu;            % time intervals of storing a snapshot
T_max = 120.0 * nu;             % maximum T for recurrent flow analysis
epsilon = 1e-6;                 % perturbation amplitude for computing J

%% initial condition
[x,~] = domain(N);              % construct the spatial domain
u0 = sin(2.0*pi*x) / sqrt(nu);  % initial condition in physical state
v0 = field2vector(u0,N);        % initial state vector

%% time integrate the KSE
[V,t_vec] = KSE_integrate(v0,t_study,dt,dt_store,nu,N);

%% computing UPOs
% to be completed


