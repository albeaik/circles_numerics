classdef DiscreteTimeEvolvingMesh < handle
    %DISCRETETIMEEVOLVINGMESH Summary of this class goes here
    %   Detailed explanation goes here
    %-- now assumes a triangulation.. can generalize to other types? 
    
    properties
        DT;
        density;
        time;
        time_step_number;
        
        x_domain_lim;
        v_domain_lim;
        x_resolution;
        v_resolution;
        
        time_history;
        DT_history;
        density_history;
        
        quality_indicators_history;
    end
    
    methods
        function InitializeMeshDefault(obj)
            %load pre specified initial conditions
            pre_specified_initial_conditions
            
            %set initial time
            time = 0;

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

            %===========| construct object
            obj.DT = DT;
            obj.density = density;
            obj.time = time;
            obj.time_step_number = 1;
            
            obj.x_domain_lim = x_domain_lim;
            obj.v_domain_lim = v_domain_lim;
            obj.x_resolution = x_resolution;
            obj.v_resolution = v_resolution;
            
            obj.time_history{obj.time_step_number} = obj.time;
            obj.DT_history{obj.time_step_number} = obj.DT;
            obj.density_history{obj.time_step_number} = obj.density;
            
            obj.quality_indicators_history = [0, 0, 0, 0, 0];
        end
        
        function EvolveMesh(obj, dt, newMesh, newMeshValues)
            obj.time = obj.time + dt;
            obj.time_step_number = obj.time_step_number + 1;
            obj.DT = newMesh;
            obj.DT_history{obj.time_step_number} = obj.DT;
            obj.density = newMeshValues;
            obj.density_history{obj.time_step_number} = obj.density;
            
            %~~~~~
            obj.MeshMaintanance();
        end
        
        function MeshMaintanance(obj)
            DT_tau = obj.DT;
            density_tau = obj.density;
            
            %~~~~~~~~ | mesh/discretization maintanance
            % (need to develop a strategy here)!! 1. evaluate, 2. adapt (in t, x, v), 3. propegate geneges or reset if necessary
            % adaptation strategy!! 1. resampling, 2. retriangulation 3. interpolation
            % for now: implement 
            %   * mesh quality and degeneracy indicators: 1. triangle areas, 2. triangle edge length, 3. triangle overlap (not implemented yet)
            %   * solution validation indicators: 1. conservation validation

            [ DT_tau_trig_areas ] = GetDelaunayAreas(DT_tau);                                           % triangle areas
            [ DT_tau_triangle_edge_lengths ] = GetDelaunayTriangleEdgeLengths( DT_tau );                    % edge lengths
            [ DT_tau_total_surface_integral ] = GetTotalDensitySurfaceIntegral(DT_tau, density_tau);    % validate conservation

            obj.quality_indicators_history(end+1, :) = [min(DT_tau_trig_areas), max(DT_tau_trig_areas), min(DT_tau_triangle_edge_lengths), max(DT_tau_triangle_edge_lengths), DT_tau_total_surface_integral];
        end
        
        function VisualizeMesh(obj)
            %~~~~~~~~ | draw --> shouldn't create side effects after this line
            visualize_trig_trig( obj.DT_history{end}, obj.density_history{end} );
            title(['Time - ', num2str(obj.time)]);
            xlim(obj.x_domain_lim*4)
            ylim(obj.v_domain_lim)
            %caxis([0, max(max(cell2mat(density_history)))])
            caxis([0, 10]);
            colorpalette = colormap('jet');
            set(gca,'Color', colorpalette(1, :));   %assuming xaxis lower limit is zero, and density of areas outside of meshed surface is also zero
        end
    end
    
end

