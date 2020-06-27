close all
clear all
clc

%load pre specified initial conditions
pre_specified_initial_conditions;
Q0 = Q0_2; %Q0_1;   %density support function

%initialize solution object: set initial domain and triangulation resolution and construct solution object
x_domain_lim = [0, 20];
v_domain_lim = [0, 20];
x_resolution = 40;
v_resolution = 40;

q_solution_mesh = DiscreteTimeEvolvingMesh(Q0, x_domain_lim, v_domain_lim, x_resolution, v_resolution);

%initialize coupled ODE object
initial_autonomous_state = [7, 7; 9, 7];  %initial state of the autonomous cars.. each row is [position, speed]
autonomous_control_u = [-2; -2];           %control for each autonomous car.. each row is [u_car_i].. assuming its constant over time for now
autonomous_car_coupling = UserDefinedCoupling(initial_autonomous_state, autonomous_control_u);

%define PDE model
PDEModel.G = @G;
PDEModel.H = @H_IDM;

%other parameters
do_realtimedraw = true;
dt = 0.01;       %t step size: seconds
T = 1;         %end time: seconds

%solve!
[solution] = characteristic_solver(q_solution_mesh, dt, T, PDEModel, autonomous_car_coupling, do_realtimedraw);

%save gif
solution.SaveAnimation(autonomous_car_coupling, '../plots/june_26/sample_sim.gif');