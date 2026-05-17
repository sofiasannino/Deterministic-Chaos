%% --- Move to this file's folder ---
cd(fileparts(mfilename('fullpath')));

%% --- Clear workspace ---
clc;
clear;
addpath("../functions/");

%% --- Output folder for plots ---
params.plot.folder = "../plots";

if ~exist(params.plot.folder, "dir")
    mkdir(params.plot.folder);
end

%% inputs
N = 128;                       % spatial resolution
epsilon = 1e-6;                % perturbation amplitude for computing J

% continuation parameters
m = 20;                        % number of continuation steps
%TBC
nu0 = 8e-4;                      % initial hyperviscosity parameter - to be completed 
nu_final = 3 * nu0 / 4;        % final hyperviscosity parameter

% initialize continuation vectors / matrices
nu_vec = linspace(nu0, nu_final, m + 1);    % vector of hyperviscosity values
v_vec = zeros(round(N / 3), m + 1);         % matrix of converged states
T_vec = zeros(1, m + 1);                    % vector of converged periods
dt_vec = 0.1 * nu_vec;                      % vector of time-steps for time integration
lambda_vec = zeros(1, m + 1);               % vector of max. Floquet multiplier

flag_vec = zeros(1, m + 1);                 % vector of convergence flags

%% converge and plot the limit cycle

%STEP1 of section 4.1:
%1)repeat section 3.1; 2)call search4PO and verify
%convergence; 3)plot in space-time and (E,P,D),

%1)
% initial condition at nu0
[x,~] = domain(N);
u_init = sin(2*pi*x) / sqrt(nu0);
v_init = field2vector(u_init,N);

% transient integration to reach the limit cycle
T_trans = 1000 * nu0;
T_int   = 250  * nu0;
dt      = dt_vec(1);


[v_trans,~] = KSE_integrate(v_init,T_trans,dt,0,nu0,N);

% store a short trajectory on the limit cycle
[v_traj,t_traj] = KSE_integrate(v_trans,T_int,dt,dt,nu0,N);

%2)
% use the last state as initial guess for the periodic point
vg0 = v_traj(:,end);

% estimate the period from consecutive maxima of the energy
[E_traj,~,~] = projection(v_traj,N,nu0);

E_traj = E_traj(:);
t_traj = t_traj(:);

imax = find( ...
    E_traj(2:end-1) > E_traj(1:end-2) & ...
    E_traj(2:end-1) > E_traj(3:end) ...
) + 1;

Tg0 = t_traj(imax(end)) - t_traj(imax(end-1));

% converge the periodic orbit
[v0,T0,flag] = search4PO(vg0,Tg0,dt,nu0,N);

% save converged solution
v_vec(:,1) = v0;
T_vec(1) = T0;
flag_vec(1) = flag; % populating flag vector for comparison at different nu

fprintf('\n Search4PO concluded with flag = %d\n',flag);
fprintf('Analysed nu = %.4e, estimated period T = %.6e  \n',nu0,T0);

% verify by integrating the converged state for one period
[v_orbit,t_orbit] = KSE_integrate(v0,T0,dt,dt,nu0,N);

%3)
% space-time plot
u_orbit = zeros(N,length(t_orbit));

for j = 1:length(t_orbit)
    u_orbit(:,j) = vector2field(v_orbit(:,j),N);
end

figure;
contourf(x,t_orbit,u_orbit.',40,'LineColor','none');
xlabel('$x$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
title(sprintf('Limit cycle at $\\nu = %.1e$',nu0),'Interpreter','latex');
colorbar;

% 3D projection in (E,P,D)
[E_orbit,P_orbit,D_orbit] = projection(v_orbit,N,nu0);

figure;
plot3(E_orbit,P_orbit,D_orbit,'LineWidth',1.5);
xlabel('$E$','Interpreter','latex');
ylabel('$P$','Interpreter','latex');
zlabel('$D$','Interpreter','latex');
grid on;
title('Limit cycle in $(E,P,D)$','Interpreter','latex');

%% first continuation step
%STEP2 of section 4.1: recompute the periodic orbit for m=20 values of nu
%going from nu0 to 3/4 nu0. use different initial guesses:
%- For i = 1, use the converged solution from ν = ν0 as the initial guess
%– For i = 2, ..., m, construct the guess state v(i)g and T(i)g by linear
%extrapolation from the two most recently converged state vectors and periods.

% first continuation guess from the converged limit cycle at nu0
i = 2;

nu = nu_vec(i);
dt = dt_vec(i);

vg = v_vec(:,i-1);
Tg = T_vec(i-1);

% converge the periodic orbit at the new value of nu
[v_new,T_new,flag] = search4PO(vg,Tg,dt,nu,N);

% save converged solution
v_vec(:,i) = v_new;
T_vec(i) = T_new;
flag_vec(i) = flag;

fprintf('\n Continuation step %d of %d , concluded with flag = %d\n',i-1,m,flag);
fprintf('Analysed nu = %.4e, estimated period T = %.6e  \n',nu,T_new);

%% subsequent continuation steps

i = 3;

while i <= m + 1

    nu = nu_vec(i);
    dt = dt_vec(i);

    % linear extrapolation from the two previous converged solutions,
    % assuming uniform spacing in nu
    vg = 2 * v_vec(:,i-1) - v_vec(:,i-2);
    Tg = 2 * T_vec(i-1)   - T_vec(i-2);

    % converge the periodic orbit at the new value of nu
    [v_new,T_new,flag] = search4PO(vg,Tg,dt,nu,N);

    % save converged solution
    v_vec(:,i) = v_new;
    T_vec(i) = T_new;
    flag_vec(i) = flag;

    fprintf('\n Continuation step %d/%d concluded with flag = %d\n',i-1,m,flag);
    fprintf('Analysed nu = %.4e, estimated period T = %.6e  \n',nu,T_new);

    i = i + 1;

end

% plot period as a function of nu
figure;
hold on;

% identify successful and unsuccessful convergence
conv_idx = flag_vec > 0;
fail_idx = flag_vec <= 0;

% plot continuation curve
plot(nu_vec,T_vec,'-','LineWidth',1.5);

% plot converged points
plot(nu_vec(conv_idx),T_vec(conv_idx),'o', ...
    'MarkerSize',6,'LineWidth',1.5);

% plot non-converged points
plot(nu_vec(fail_idx),T_vec(fail_idx),'x', ...
    'MarkerSize',8,'LineWidth',1.5);

xlabel('$\nu$','Interpreter','latex');
ylabel('$T$','Interpreter','latex');
title('Continuation of the limit cycle period','Interpreter','latex');
legend('continuation branch','converged','not converged', ...
    'Location','best','Interpreter','latex');
grid on;

% print values required in the report table
selected_idx = [1 6 11 16 21];

fprintf('\n ---- Values for the report table ----:\n');
fprintf('      nu              T              flag\n');

for j = selected_idx
    fprintf(' %.6e    %.6e      %d\n',nu_vec(j),T_vec(j),flag_vec(j));
end




%% Floquet exponents

%STEP3 of section 4.1: 1)repeat the computation of step 3.2 for each nu and
%integrate for exactly one trajectory; 2)plot Floquet multipliers as
%function of nu e report the relevant ones in the table

for i = 1:(m + 1)

    nu = nu_vec(i);
    dt = dt_vec(i);
    T  = T_vec(i);
    v0 = v_vec(:,i);

    % compute the Jacobian of the flow map over one period
    J = Jacobian(v0,T,epsilon,dt,nu,N);

    % Floquet multipliers are the eigenvalues of the Jacobian
    lambda = eig(J);

    % remove the trivial multiplier lambda = 1
    lambda_tol = 1e-6;
    idx_trivial = find(abs(lambda - 1) < lambda_tol);
    
    if isempty(idx_trivial)
    
        % if lambda = 1 is not exact because of numerical error,
        % remove the closest multiplier to 1
        [~,idx_trivial] = min(abs(lambda - 1));
    
    end

    lambda(idx_trivial) = [];

    % store the largest non-trivial multiplier modulus
    lambda_vec(i) = max(abs(lambda));

    fprintf('\n Floquet calculation %d/%d completed\n',i,m+1);
    fprintf('Analysed nu = %.4e, |lambda*| = %.6e  \n',nu,lambda_vec(i));

end

% plot the largest non-trivial Floquet multiplier as a function of nu
figure;
plot(nu_vec,lambda_vec,'o-','LineWidth',1.5);
hold on;
yline(1,'--','LineWidth',1.2);
xlabel('$\nu$','Interpreter','latex');
ylabel('$|\lambda^\ast|$','Interpreter','latex');
title('Floquet stability of the continued limit cycle','Interpreter','latex');
grid on;

% estimate the crossing of |lambda*| = 1
idx_cross = find((lambda_vec(1:end-1)-1).*(lambda_vec(2:end)-1) <= 0,1,'first');

if ~isempty(idx_cross)

    nu_cross = nu_vec(idx_cross) + ...
        (1 - lambda_vec(idx_cross)) * ...
        (nu_vec(idx_cross+1) - nu_vec(idx_cross)) / ...
        (lambda_vec(idx_cross+1) - lambda_vec(idx_cross));

    fprintf('\n Estimated crossing of |lambda*| = 1:\n');
    fprintf('nu_cross = %.6e\n',nu_cross);

end

% print values required in the report table
selected_idx = [1 6 11 16 21];

fprintf('\n ----- Values for the report table: -----\n');
fprintf('      nu              |lambda*|\n');

for j = selected_idx
    fprintf(' %.6e    %.6e\n',nu_vec(j),lambda_vec(j));
end
