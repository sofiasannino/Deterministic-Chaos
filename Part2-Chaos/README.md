# EPFL ME467: Turbulence
### Time-marching the KSE
The following lines show how to utilize the function `KSE_integrate` to advance the KSE in time.
```
close all; clear; clc

addpath('../functions/')
nu = 1e-2;
N = 128;
T = 750 * nu;
dt = 0.1 * nu;

[x,~] = domain(N);
u0 = sin(2.0*pi*x) / sqrt(nu);

v0 = field2vector(u0,N);
[vT,~] = KSE_integrate(v0,T,dt,0,nu,N);

uT = vector2field(vT,N);

figure
  plot(x,uT,'LineWidth',2)
  hold on; grid on
  plot(x,u0,'LineWidth',2)
  xlabel('x'); ylabel('u')
  legend('u(T)','u(0)')
```
This program can be broken down into the following steps:
- We first tell MATLAB to look for functions in the relative folder path `../functions/`.
- Then, we define the hyperviscosity parameter `nu = 1e-2`, number of grid points `N = 128`, integration time `T = 750 * nu` and time step size `dt = 0.1 * nu`.
- To define the initial condition, we first construct the grid points using `[x,~] = domain(N)`. Then, the initial condition is defined using the sine function `u0 = sin(2.0*pi*x) / sqrt(nu)`. Note that this initial condition is consistent with the imposed center symmetry.
- We construct the state vector associated with `u0` by calling `v0 = field2vector(u0,N)`. Finally, we call `[vT,~] = KSE_integrate(v0,T,dt,0,nu,N)` where `vT` is the state vector of the KSE obtained by advancing `v0` for time `T`.
- Since `vT` is a state vector and not the vector of grid values, we call `uT = vector2field(vT,N)` to transform `vT` to its corresponding physical field `uT` and, then, plot `u0` and `uT` in one figure.

 **Note 1**: By setting the 4th argument of `KSE_integrate` to `0`, no intermediate snapshot is stored and `vT` is a column vector. In order to store intermediate snapshots, you can call `[vt,t] = KSE_integrate(v0,T,dt,dt_store,nu,N)`. As a result, every `dt_store` the snapshot is appended to `vt` as a new column. The vector `t` contains the time instance associated with each column of `vt`.

### 3D projection of the state space
The following lines show how to use the function `projection` to visualize a 3D projection of the chaotic trajectory from the previous example.
```
close all; clear; clc

addpath('../functions/')
nu = 6.5e-4;
N = 128;
T = 750 * nu;
dt = 0.1 * nu;
dt_store = 1 * nu;

[x,~] = domain(N);
u0 = sin(2*pi*x) / sqrt(nu);

v0 = field2vector(u0,N);
[vv,tt] = KSE_integrate(v0,T,dt,dt_store,nu,N);

[E,P,D] = projection(vv,N,nu);

figure
  plot3(E,P,D,'Color',[0,0,0,0.25])
  xlabel('E'); ylabel('P'); zlabel('D'); grid on
```
Here `[0,0,0,0.25]` specifies red, green, blue and the opacity of the line. We have chosen a semi-transparent line so that regions of the state space that are visited more often appear darker than the regions that are rarely visited.

### Search for a periodic orbit 
The following program is an example of how to utilize the function `search4PO` to compute a periodic orbit from a guess. Here, the same sine function as in the previous examples is used, and a period of `T_guess=50 * nu` is chosen for demonstration.
```
close all; clear; clc

addpath('../functions/')
nu = 6.5e-4;
N = 128;
dt = 0.1 * nu;

[x,~] = domain(N);
u0 = sin(2*pi*x) / sqrt(nu);

v_guess = field2vector(u0,N);
T_guess = 50 * nu;

[v_best,T_best,flag] = search4PO(v_guess,T_guess,dt,nu,N); 
```
This crude guess does not converge, and the following will be displayed:
```
                                         Norm of      First-order   Trust-region
 Iteration  Func-count     f(x)          step         optimality    radius
     0         45     3.66815e+07                      5.72e+08               1
     1         90     3.53354e+07              1       2.23e+09               1
   ...        ...             ...            ...            ...             ...
    75       3244     2.44136e+07      0.0610352       1.56e+05           0.061

Solver stopped prematurely.

fsolve stopped because it exceeded the iteration limit,
options.MaxIterations = 7.500000e+01.
```
In a successful search, the search flag takes the value of 1. Here, however, `flag = 0`. If the search is successful, `v_best` is the converged state vector, and `T_best` is the converged period of the periodic orbit. In that case, if you advance `v_best` as the initial condition for a time interval of `T_best` (see the first example above), the trajectory closes on itself. Note that `v_guess` and `v_best` are both state vectors.
