close all; clear; clc

%% inputs
addpath('../functions/')       % folder containing functions
nu = 1e-2;                     % hyperviscosity parameter
N = 128;                       % spatial resolution
T_trans = 1000 * nu;           % transient time period
T_int = 250 * nu;              % analysis time period
dt = 0.1 * nu;                 % time step size for time integration
dt_store = 0.1 * nu;           % time intervals of storing a snapshot
epsilon = 1e-2;                % relative amplitude of perturbation

%% Step 1: initial condition
[x,~] = domain(N);             % construct the spatial domain
u0 = sin(2.0*pi*x) / sqrt(nu); % initial condition in physical state
v0 = field2vector(u0,N);       % initial state vector

%% Step 2: transient time integration
[v1000,~] = KSE_integrate(v0,T_trans,dt,0,nu,N);

%% Step 3: perturbing the state vector
r = zeros(size(v1000));        % the perturbation vector memory allocation
for k = 1:length(r)
    %%% to be completed
end

%%% to be completed
