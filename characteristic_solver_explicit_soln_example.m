close all
clear all
clc

%load pre specified initial conditions
pre_specified_initial_conditions;
Q0 = Q0_5; Q0_2; %Q0_1;   %density support function

%Q0 = @(x,v) (x >= 8) .* (x <= 12);


resolution_set = [20, 40, 60];
dt_set = [0.1, 0.05, 0.01];

%initialize solution object: set initial domain and triangulation resolution and construct solution object

solutions = [];
index = 1;
for res = resolution_set
    for dt = dt_set
        %----
        initialization_option = 'oracle-function';
        initialization_parameters.Q0 = @(x,v) (x >= -0.5) & (x <= 0.5) & (v >= -0.5) & (v <= 0.5);
        %initialization_parameters.mesh_domain_limits = [-0.5, 0.5; -0.5, 0.5];
        %initialization_parameters.mesh_domain_limits = [-1, 1; -1, 1];
        initialization_parameters.mesh_domain_limits = [-0.6, 0.6; -0.6, 0.6];
        initialization_parameters.mesh_resolution = res*[1; 1];
        q_solution_mesh = DiscreteTimeEvolvingMesh(initialization_option, initialization_parameters);

        %initialize coupled ODE object ~ ignored
        activate_coupling = false;
        initial_autonomous_state = [7, 7; 9, 7];  %initial state of the autonomous cars.. each row is [position, speed]
        autonomous_control_u = [-2; -2];           %control for each autonomous car.. each row is [u_car_i].. assuming its constant over time for now
        autonomous_car_coupling = UserDefinedCoupling(activate_coupling, initial_autonomous_state, autonomous_control_u);

        %define PDE model
        PDEModel.G = @H_explicit_soln_example;
        PDEModel.H = @H_explicit_soln_example;

        %other parameters
        do_realtimedraw = false;
        %dt = 0.01;       %t step size: seconds
        T = 4;         %end time: seconds

        %solve!
        tic;
        solutions{index}.q = characteristic_solver(q_solution_mesh, dt, T, PDEModel, autonomous_car_coupling, do_realtimedraw);
        solutions{index}.computational_time = toc;
        solutions{index}.coupling = autonomous_car_coupling;
        solutions{index}.res = res;
        solutions{index}.dt = dt;
        index = index + 1;

        %visualizations
        %solution.VisualizeAnimation(autonomous_car_coupling, 2)
        %solution.SaveAnimation(autonomous_car_coupling, '../plots/june_26/sample_sim.gif');
    end
end

%comparison with explicit solution
xi_explicit_soln = @(t) nansum([(t<1), (t>1 & t<=2), (t>2)] .* [(1/(2-t)-1/2), (t-1/2), (2.5-1/(t-1))]);
q_explicit_soln = @(t, x1, x2) (x1 >= xi_explicit_soln(t)-1/2) & (x1 <= xi_explicit_soln(t)+1/2) & (x2 >= xi_explicit_soln(t)-1/2) & (x2 <= xi_explicit_soln(t)+1/2);

soln_dif_t = [];
for soln_number = 1:size(solutions, 2)
    solution_package = solutions{soln_number};
    solution = solution_package.q;
    soln_dif_t{soln_number} = zeros([solution.time_step_number, 1]);
    for t_index = 1:solution.time_step_number
        t = solution.time_history{t_index};
        numerical_density = solution.density_history{t_index};

        num_DT = solution.DT_history{t_index};

        centroids = solution.GetMeshCentroids(num_DT);
        explicit_density = q_explicit_soln(t, centroids(:, 1), centroids(:, 2));

        areas = solution.GetMeshAreas(num_DT);

        density_error = numerical_density - explicit_density;
        soln_dif_t{soln_number}(t_index) = sqrt(sum((density_error.*areas).^2));

        %----
        %solution.visualize_triangulation( num_DT, density_error );
        %title(['Time - ', num2str(t)]);
        %xlim([-1, 5])
        %ylim([-1, 5])
        %%caxis([0, max(max(cell2mat(density_history)))])
        %caxis([0, 2]);
        %colorpalette = colormap('jet');
        %%set(gca,'Color', colorpalette(1, :));   %assuming xaxis lower limit is zero, and density of areas outside of meshed surface is also zero
        %drawnow
    end
end

figure
for soln_number = 1:size(solutions, 2)
    subplot(size(resolution_set, 2), size(dt_set, 2), soln_number)
    
    solution_package = solutions{soln_number};
    solution = solution_package.q;
    mesh_res = solution_package.res;
    dt = solution_package.dt;
    error = soln_dif_t{soln_number};
    
    plot(cell2mat(solution.time_history), error, 'linewidth', 2)
    xlim([0, T])
    ylim([0, 0.02])
    %xlabel('t')
    %ylabel('Density error [L^2 mass]')
    %title(['res - ', num2str(mesh_res), ' - dt - ', num2str(dt)])
    set(gca,'FontSize',15);
    grid on
    
    if(soln_number <= size(dt_set, 2))
        title(['dt = ', num2str(dt)])
    end
    
    if(mod(soln_number-1, size(dt_set, 2))+1 == 1)
        ylabel(['Mesh resolution  = ', num2str(mesh_res)])
    end
    
    if(soln_number > size(resolution_set, 2)*size(dt_set, 2)-size(dt_set, 2))
        xlabel('t')
    end
end

runtimes = zeros([size(resolution_set, 2), size(dt_set, 2)]);
for soln_number = 1:size(solutions, 2)
    solution_package = solutions{soln_number};
    solution = solution_package.q;
    mesh_res = solution_package.res;
    dt = solution_package.dt;
    runtime = solution_package.computational_time;
    
    
    runtimes(ceil(soln_number/size(resolution_set, 2)), mod(soln_number-1, size(dt_set, 2))+1) = runtime;
end


return

solution = solutions{index}.q;

%compute l2 mass of error --------------------------
soln_dif_t = zeros([solution.time_step_number, 1]);
for t_index = 1:solution.time_step_number
    t = solution.time_history{t_index};
    numerical_density = solution.density_history{t_index};
    
    num_DT = solution.DT_history{t_index};
    
    centroids = solution.GetMeshCentroids(num_DT);
    explicit_density = q_explicit_soln(t, centroids(:, 1), centroids(:, 2));
    
    areas = solution.GetMeshAreas(num_DT);
    
    density_error = numerical_density - explicit_density;
    soln_dif_t(t_index) = sqrt(sum((density_error.*areas).^2));
    
    %----
    %solution.visualize_triangulation( num_DT, density_error );
    %title(['Time - ', num2str(t)]);
    %xlim([-1, 5])
    %ylim([-1, 5])
    %%caxis([0, max(max(cell2mat(density_history)))])
    %caxis([0, 2]);
    %colorpalette = colormap('jet');
    %%set(gca,'Color', colorpalette(1, :));   %assuming xaxis lower limit is zero, and density of areas outside of meshed surface is also zero
    %drawnow
end

plot(cell2mat(solution.time_history), soln_dif_t, 'linewidth', 2)
xlim([0, T])
xlabel('t')
ylabel('Density error [L^2 mass]')
title('Numerical solution compared to explicit solution')
set(gca,'FontSize',15);
grid on

%animate explicit and numerical solutions --------------------------
for t_index = 1:solution.time_step_number
    t = solution.time_history{t_index};
    clf
    
    numerical_density = solution.density_history{t_index};
    num_DT = solution.DT_history{t_index};
    
    explicit_soln_boudnary = xi_explicit_soln(t);
    %rectangle('Position', [explicit_soln_boudnary-0.5, explicit_soln_boudnary-0.5 explicit_soln_boudnary+0.5, explicit_soln_boudnary+0.5], 'EdgeColor','r');
    plot3( [explicit_soln_boudnary-0.5 explicit_soln_boudnary+0.5 explicit_soln_boudnary+0.5 explicit_soln_boudnary-0.5 explicit_soln_boudnary-0.5], [explicit_soln_boudnary-0.5 explicit_soln_boudnary-0.5 explicit_soln_boudnary+0.5 explicit_soln_boudnary+0.5 explicit_soln_boudnary-0.5], 10*[1 1 1 1 1], 'Color','r')
    hold on
    
    solution.visualize_triangulation( num_DT, numerical_density );
    title(['Time - ', num2str(t)]);
    xlim([-1, 3])
    ylim([-1, 3])
    %caxis([0, max(max(cell2mat(density_history)))])
    caxis([0, 2]);
    colorpalette = colormap('jet');
    %set(gca,'Color', colorpalette(1, :));   %assuming xaxis lower limit is zero, and density of areas outside of meshed surface is also zero
    drawnow
end


%animate explicit and numerical solutions --------------------------
fighandle = figure
plt_time_instances = [0, 1, 2, 3, 4]./dt+1;
for index = 1:size(plt_time_instances, 2)
    subplot(size(plt_time_instances, 2), 1, index)
    
    t_index = plt_time_instances(index);
    t = solution.time_history{t_index};
    
    numerical_density = solution.density_history{t_index};
    num_DT = solution.DT_history{t_index};
    
    explicit_soln_boudnary = xi_explicit_soln(t);
    %rectangle('Position', [explicit_soln_boudnary-0.5, explicit_soln_boudnary-0.5 explicit_soln_boudnary+0.5, explicit_soln_boudnary+0.5], 'EdgeColor','r');
    plot3( [explicit_soln_boudnary-0.5 explicit_soln_boudnary+0.5 explicit_soln_boudnary+0.5 explicit_soln_boudnary-0.5 explicit_soln_boudnary-0.5], [explicit_soln_boudnary-0.5 explicit_soln_boudnary-0.5 explicit_soln_boudnary+0.5 explicit_soln_boudnary+0.5 explicit_soln_boudnary-0.5], 10*[1 1 1 1 1], 'Color','r', 'linewidth', 1)
    hold on
    
    solution.visualize_triangulation( num_DT, numerical_density );
    title(['Time - ', num2str(t)]);
    xlim([-1, 3])
    ylim([-1, 3])
    %caxis([0, max(max(cell2mat(density_history)))])
    caxis([0, 1]);
    colorpalette = colormap(flipud(gray));
    set(gca,'Color', colorpalette(1, :));   %assuming xaxis lower limit is zero, and density of areas outside of meshed surface is also zero
    
    grid on
end

set(gcf, 'Position',  [560    68   231   880])


%plot error between explicit and numerical solutions ---------------------
fighandle = figure
plt_time_instances = [0, 1, 2, 3, 4]./dt+1;
for index = 1:size(plt_time_instances, 2)
    subplot(size(plt_time_instances, 2), 1, index)
    
    t_index = plt_time_instances(index);
    t = solution.time_history{t_index};
    numerical_density = solution.density_history{t_index};
    
    num_DT = solution.DT_history{t_index};
    
    centroids = solution.GetMeshCentroids(num_DT);
    explicit_density = q_explicit_soln(t, centroids(:, 1), centroids(:, 2));
    
    areas = solution.GetMeshAreas(num_DT);
    
    density_error = numerical_density - explicit_density;
    soln_dif_t(t_index) = sqrt(sum((density_error.*areas).^2));
    
    %----
    solution.visualize_triangulation( num_DT, abs(density_error) );
    title(['Time - ', num2str(t)]);
    xlim([-1, 3])
    ylim([-1, 3])
    %caxis([0, max(max(cell2mat(density_history)))])
    caxis([0, 1]);
    colorpalette = colormap(flipud(gray));
    set(gca,'Color', colorpalette(1, :));   %assuming xaxis lower limit is zero, and density of areas outside of meshed surface is also zero
    drawnow
end

set(gcf, 'Position',  [560    68   231   880])