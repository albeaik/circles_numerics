close all
clear all
clc


do_realtimedraw = true;
dt = 0.01;       %t step size: seconds
T = 1.5;         %end time: seconds

autonomous_car_coupling = UserDefinedCoupling;
autonomous_car_coupling.InitializeDefault();

solution_mesh = DiscreteTimeEvolvingMesh;
solution_mesh.InitializeMeshDefault();

[solution] = characteristic_solver(solution_mesh, dt, T, @G, @H_IDM, autonomous_car_coupling, do_realtimedraw);