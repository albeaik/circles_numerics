close all
clear all
clc

%load pre specified initial conditions
pre_specified_initial_conditions;
Q0 = Q0_2; %Q0_1;   %density support function

%initialize solution object: set initial domain and triangulation resolution and construct solution object

%----
% initialization_option = 'oracle-function';
% initialization_parameters.Q0 = Q0;
% initialization_parameters.mesh_domain_limits = [0, 20; 0, 20];
% initialization_parameters.mesh_resolution = [40; 40];
% q_solution_mesh = DiscreteTimeEvolvingMesh(initialization_option, initialization_parameters);

%----
load('sample_mesh_and_density')
%pdecirc(10,10, 5)
%pderect([0 20 0 20])
initialization_option = 'pre-defined-initialization';
initialization_parameters.DT = DT;
initialization_parameters.density = density;
q_solution_mesh = DiscreteTimeEvolvingMesh(initialization_option, initialization_parameters);

%initialize coupled ODE object
activate_coupling = true;
initial_autonomous_state = [7, 7; 9, 7];  %initial state of the autonomous cars.. each row is [position, speed]
autonomous_control_u = [-2; -2];           %control for each autonomous car.. each row is [u_car_i].. assuming its constant over time for now
autonomous_car_coupling = UserDefinedCoupling(activate_coupling, initial_autonomous_state, autonomous_control_u);

%define PDE model
PDEModel.G = @G;
PDEModel.H = @H_IDM;

%other parameters
do_realtimedraw = true;
dt = 0.01;       %t step size: seconds
T = 1;         %end time: seconds

%solve!
[solution] = characteristic_solver(q_solution_mesh, dt, T, PDEModel, autonomous_car_coupling, do_realtimedraw);

%visualizations
solution.VisualizeAnimation(autonomous_car_coupling, 2)
%solution.SaveAnimation(autonomous_car_coupling, '../plots/june_26/sample_sim.gif');
%VisualizeMeshEvolvingInTime(solution, autonomous_car_coupling ,[0, 1])
%solution.VisualizeMultipleSteps(autonomous_car_coupling, [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])