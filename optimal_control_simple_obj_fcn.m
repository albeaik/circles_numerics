function [ obj_value ] = optimal_control_simple_obj_fcn( u, nAC, nT )
%OPTIMAL_CONTROL_SIMPLE_OBJ_FCN Summary of this function goes here
%   Detailed explanation goes here

    %obj_value = sum(u.^2);
    

    %load pre specified initial conditions
    pre_specified_initial_conditions;
    Q0 = Q0_2; %Q0_1;   %density support function

    %initialize solution object: set initial domain and triangulation resolution and construct solution object
    mesh_domain_limits = [0, 20; 0, 20];
    mesh_resolution = [40; 40];
    q_solution_mesh = DiscreteTimeEvolvingMesh(Q0, mesh_domain_limits, mesh_resolution);

    %initialize coupled ODE object
    activate_coupling = true;
    initial_autonomous_state = [7, 7; 9, 7];    %initial state of the autonomous cars.. each row is [position, speed]
    autonomous_control_u = reshape(u, nAC, nT);                   %control for each autonomous car.. each row is [u_car_i].. assuming its constant over time for now
    autonomous_car_coupling = UserDefinedCoupling(activate_coupling, initial_autonomous_state, autonomous_control_u);

    %define PDE model
    PDEModel.G = @G;
    PDEModel.H = @H_IDM;

    %other parameters
    do_realtimedraw = false;
    dt = 0.01;       %t step size: seconds
    T = dt*nT;         %end time: seconds

    %solve!
    [solution] = characteristic_solver(q_solution_mesh, dt, T, PDEModel, autonomous_car_coupling, do_realtimedraw);

    obj_value = 0;
    for i = 1:solution.time_step_number
        centroids = DiscreteTimeEvolvingMesh.GetMeshCentroids( solution.DT_history{i} );
        obj_value = obj_value + (solution.density_history{i}' * centroids(:, 2)) ./ size(centroids, 1);
    end
    
    obj_value = obj_value ./ solution.time_step_number
end

