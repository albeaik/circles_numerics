function [solution] = characteristic_solver(G_fcn, H_fcn)
    %load pre specified initial conditions
    pre_specified_initial_conditions

    %identify initial domain
    x_domain_lim = [0, 20];
    v_domain_lim = [0, 20];
    x_resolution = 40;
    v_resolution = 40;

    %n = 500;

    %--- std basic meshing
    % init_domain_points = [x_domain_lim(1), v_domain_lim(1); x_domain_lim(2), v_domain_lim(1); x_domain_lim(2), v_domain_lim(2); x_domain_lim(1), v_domain_lim(2); x_domain_lim(1), v_domain_lim(1)];   % (x, v)
    % %init_domain_points = [init_domain_points; rand(n, 1)*x_domain_lim(2), rand(n, 1)*v_domain_lim(2)];
    % 
    % [p,t,e] = pmesh(init_domain_points, 50, 6);
    % %e=boundary_nodes(t);
    % 
    % DT.Points = p;
    % DT.ConnectivityList = t;

    %--- matlab native meshing
    x=linspace(x_domain_lim(1), x_domain_lim(2), x_resolution);
    v=linspace(v_domain_lim(1), v_domain_lim(2), v_resolution);
    [x_grid, v_grid]=meshgrid(x, v);

    init_domain_points = [x_grid(:), v_grid(:)];

    DT = delaunayTriangulation(init_domain_points);
    p = DT.Points;
    t = DT.ConnectivityList;


    %--- generate density profile
    %Q0 = Q0_1;   %density support function
    Q0 = Q0_2;

    tmp_trig_centers = DT.incenter;
    %tmp_trig_centers = GetDelaunayCentroids( DT );
    density = zeros([size(DT.ConnectivityList, 1), 1]);
    density(Q0(tmp_trig_centers(:, 1), tmp_trig_centers(:, 2)) == 1) = 1;

    %tmp_trig_centers = DT.incenter;
    %tmp_trig_centers = GetDelaunayCentroids( DT );
    %triplot(DT)
    %trisurf(DT.ConnectivityList, DT.Points(:, 1), DT.Points(:, 2), density)
    %plot3(tmp_trig_centers(:, 1), tmp_trig_centers(:, 2), density, '.', 'markersize', 15)

    %--- identify discontinuities
    % assuming discontinuities were not known apriori, but are detectable
    maxdtdensity = 0.5;
    discontinuity_edge_list = get_discontinuity_edge_list(DT, density, maxdtdensity);

    DT.Constraints = discontinuity_edge_list;   %add discontinuity_edge_list as a constraint to triangulation

    %--- solve

    time = 0;       %time: seconds
    dt = 0.01;       %t step size: seconds
    T = 3;          %end time: seconds
    remesh_cycle = 100;  %time steps
    assume_fixed_DT = 0;    %used to determain if pre-computations can be used ==> under what conditions can I assume fixed DT?! need to confirm
    mesh_subarea_minmax_limit = [0.0005, 2]; %min and max allowed triangle area before mesh refinement
    dt_minmax_limit = [0.01, 0.0001];
    refine_mesh_now = 0;

    Autonomous_init_state = [7, 7; 9, 7];  %initial state of the autonomous cars.. each row is [position, speed]
    Autonomous_control_u = [-2; -2];           %control for each autonomous car.. each row is [u_car_i].. assuming its constant over time for now

    tmp_quality_indicators_history = [0, 0, 0, 0, 0];

    DT_t = DT;
    DT_history{1} = DT_t;
    density_t = density;
    density_history{1} = density_t;
    Autonomous_state_t = Autonomous_init_state;
    Autonomous_state_history{1} = Autonomous_state_t;
    t_index = 1;
    time = 0;
    %for time = 0:dt:T                   %NOTE: inside this loop, DT_.. changes from native native to non-native but compatible data type
    while(time <= T)
    % write down iteration contract: assume DT_t, density_t set properly, never changed until end of iteration, etc

        %~~~~~~~~ | prepare for this iteration
        time
        %DT_tau = DT_t;
        %density_tau = density_t;

        [DT_t_centroids] = GetDelaunayCentroids(DT_t);  %sorted by triangle id --> can pre-compute if DT is fixed
        [DT_t_areas] = GetDelaunayAreas(DT_t);  %sorted by triangle id --> can pre-compute if DT is fixed

        %~~~~~~~~ | evolve solution in time
        % simulate the characteristic equation
        DT_tau_Points = DT_t.Points + dt*[G_fcn(DT_t, density_t, Autonomous_state_t, Autonomous_state_t, assume_fixed_DT), H_fcn(DT_t, density_t, Autonomous_state_t, DT_t.Points, assume_fixed_DT)];

        % simulate the autonomous car ODE
        Autonomous_state_tau = Autonomous_state_t + [Autonomous_state_t(:, 2)*dt, (H_IDM(DT_t, density_t, Autonomous_state_t, Autonomous_state_t, assume_fixed_DT) + Autonomous_control_u) * dt];

        %create DT_tau object ==> note matlab native DT object recomputes triangulation if points updated directly                       
        DT_tau.Points = DT_tau_Points;
        DT_tau.ConnectivityList = DT_t.ConnectivityList;

        % simulate density scaling
        [DT_t_trig_areas] = GetDelaunayAreas(DT_t);  %sorted by triangle id --> can pre-compute if DT is fixed
        [DT_tau_trig_areas] = GetDelaunayAreas(DT_tau);  %sorted by triangle id
        density_tau = (DT_t_trig_areas ./ DT_tau_trig_areas) .* density_t;

        %~~~~~~~~ | mesh/discretization maintanance
        % (need to develop a strategy here)!! 1. evaluate, 2. adapt (in t, x, v), 3. propegate geneges or reset if necessary
        % adaptation strategy!! 1. resampling, 2. retriangulation 3. interpolation
        % for now: implement 
        %   * mesh quality and degeneracy indicators: 1. triangle areas, 2. triangle edge length, 3. triangle overlap (not implemented yet)
        %   * solution validation indicators: 1. conservation validation

        [ DT_tau_trig_areas ] = GetDelaunayAreas(DT_tau);                                           % triangle areas
        [ DT_tau_triangle_edge_lengths ] = GetDelaunayTriangleEdgeLengths( DT );                    % edge lengths
        [ DT_tau_total_surface_integral ] = GetTotalDensitySurfaceIntegral(DT_tau, density_tau);    % validate conservation

        tmp_quality_indicators_history(end+1, :) = [min(DT_tau_trig_areas), max(DT_tau_trig_areas), min(DT_tau_triangle_edge_lengths), max(DT_tau_triangle_edge_lengths), DT_tau_total_surface_integral];

        %~~~~~~~~ | prepare for next iteration
        time = time + dt;
        t_index = t_index + 1;
        DT_t = DT_tau; % DT_tau_remished; % DT_tau;
        DT_history{t_index} = DT_t;
        density_t = density_tau; % density_DT_tau_remished; % density_tau;
        density_history{t_index} = density_t;
        Autonomous_state_t = Autonomous_state_tau;
        Autonomous_state_history{t_index} = Autonomous_state_t;

        %~~~~~~~~ | draw --> shouldn't create side effects after this line
        visualize_trig_trig( DT_history{end}, density_history{end} );
        title(['time - ', num2str(time)]);
        xlim(x_domain_lim*4)
        ylim(v_domain_lim)
        %caxis([0, max(max(cell2mat(density_history)))])
        caxis([0, 5]);
        colorpalette = colormap('jet');
        set(gca,'Color', colorpalette(1, :));   %assuming xaxis lower limit is zero, and density of areas outside of meshed surface is also zero

        hold on
        plot3(Autonomous_state_history{end}(:, 1), Autonomous_state_history{end}(:, 2), 10*ones(size(Autonomous_state_history{end}(:, 2))), 'r*')
        hold off

        drawnow
    end
    
    solution = density_history;
end