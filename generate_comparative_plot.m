close all
clear all
clc

%load pre specified initial conditions
pre_specified_initial_conditions;
Q0 = Q0_3; Q0_2; %Q0_1;   %density support function

%variations
mesh_resolution_set = [40, 20];
dt_set = [0.1, 0.01];

index = 1;
solutions = [];
for mesh_res = mesh_resolution_set
    for dt = dt_set
        %initialize solution object: set initial domain and triangulation resolution and construct solution object
        initialization_option = 'oracle-function';
        initialization_parameters.Q0 = Q0;
        initialization_parameters.mesh_domain_limits = [0, 20; 0, 20];
        initialization_parameters.mesh_resolution = mesh_res*[1; 1];
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
        do_realtimedraw = false;
        %dt = 0.01;       %t step size: seconds
        T = 1;         %end time: seconds

        %solve!
        tic;
        solutions{index}.q = characteristic_solver(q_solution_mesh, dt, T, PDEModel, autonomous_car_coupling, do_realtimedraw);
        solutions{index}.computational_time = toc;
        solutions{index}.coupling = autonomous_car_coupling;
        solutions{index}.res = mesh_res;
        solutions{index}.dt = dt;
        index = index + 1;

        %visualizations
        %solution.VisualizeAnimation(autonomous_car_coupling, 2)
        %solution.SaveAnimation(autonomous_car_coupling, '../plots/june_26/sample_sim.gif');
    end
end

figure
plt_times = [0, 0.5, 1];
for index = 1:size(solutions, 2)
    for t_ind = 1:size(plt_times, 2)
        subplot(size(solutions, 2), size(plt_times, 2), (index-1)*size(plt_times, 2)+t_ind)
        time = plt_times(t_ind);
        n = plt_times(t_ind)/solutions{index}.dt+1; %find(cell2mat(solutions{index}.q.time_history) == time);
        solutions{index}.q.VisualizeStep(n);
        solutions{index}.coupling.VisualizeStep(n);
        
        if(index ~= 1)
            title('')
        end
        
        if(t_ind == 1)
            ylabel(['mesh ', num2str(solutions{index}.res), ' - dt ', num2str(solutions{index}.dt)])
        end
    end
end
