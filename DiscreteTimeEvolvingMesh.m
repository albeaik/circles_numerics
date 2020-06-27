classdef DiscreteTimeEvolvingMesh < handle
    %DISCRETETIMEEVOLVINGMESH Summary of this class goes here
    %   Detailed explanation goes here
    %-- now assumes a triangulation.. can generalize to other types? 
    
    properties
        DT;
        density;
        time;
        time_step_number;
        
        mesh_domain_limits; %[x_domain_lim; v_domain_lim]   %2 by 2 matrix
        mesh_resolution;    %[x_resolution; v_resolution]   %2 by 1 vector
        
        time_history;
        DT_history;
        density_history;
        
        quality_indicators_history;
    end
    
    methods
        function obj = DiscreteTimeEvolvingMesh(Q0, mesh_domain_limits, mesh_resolution)
            %default constuctor function!
            
            %set initial time
            time = 0;
            
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
            x=linspace(mesh_domain_limits(1, 1), mesh_domain_limits(1, 2), mesh_resolution(1));
            v=linspace(mesh_domain_limits(2, 1), mesh_domain_limits(2, 2), mesh_resolution(2));
            [x_grid, v_grid]=meshgrid(x, v);

            init_domain_points = [x_grid(:), v_grid(:)];

            DT = delaunayTriangulation(init_domain_points);
            %p = DT.Points;
            %t = DT.ConnectivityList;


            %--- generate density profile
            %tmp_trig_centers = DT.incenter;
            tmp_trig_centers = DiscreteTimeEvolvingMesh.GetMeshCentroids( DT );
            density = zeros([size(DT.ConnectivityList, 1), 1]);
            density(Q0(tmp_trig_centers(:, 1), tmp_trig_centers(:, 2)) == 1) = 1;

            %tmp_trig_centers = DT.incenter;
            %tmp_trig_centers = DiscreteTimeEvolvingMesh.GetMeshCentroids( DT );
            %triplot(DT)
            %trisurf(DT.ConnectivityList, DT.Points(:, 1), DT.Points(:, 2), density)
            %plot3(tmp_trig_centers(:, 1), tmp_trig_centers(:, 2), density, '.', 'markersize', 15)

            %--- identify discontinuities
            % assuming discontinuities were not known apriori, but are detectable
            %maxdtdensity = 0.5;
            %discontinuity_edge_list = get_discontinuity_edge_list(DT, density, maxdtdensity);

            %DT.Constraints = discontinuity_edge_list;   %add discontinuity_edge_list as a constraint to triangulation

            %===========| construct object
            obj = obj@handle();
            
            obj.DT = DT;
            obj.density = density;
            obj.time = time;
            obj.time_step_number = 1;
            
            %obj.x_domain_lim = x_domain_lim;
            %obj.v_domain_lim = v_domain_lim;
            %obj.x_resolution = x_resolution;
            %obj.v_resolution = v_resolution;
            
            obj.mesh_domain_limits = mesh_domain_limits;
            obj.mesh_resolution = mesh_resolution;
            
            obj.time_history{obj.time_step_number} = obj.time;
            obj.DT_history{obj.time_step_number} = obj.DT;
            obj.density_history{obj.time_step_number} = obj.density;
            
            obj.quality_indicators_history = [0, 0, 0, 0, 0];
        end
        
        function EvolveMesh(obj, dt, newMesh, newMeshValues)
            obj.time_step_number = obj.time_step_number + 1;
            
            obj.time = obj.time + dt;
            obj.time_history{obj.time_step_number} = obj.time;
            
            obj.DT = newMesh;
            obj.DT_history{obj.time_step_number} = obj.DT;
            
            obj.density = newMeshValues;
            obj.density_history{obj.time_step_number} = obj.density;
            
            %~~~~~
            obj.MeshMaintanance();  %now only updates quality metrics.. figure out what to do here later
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

            [ DT_tau_trig_areas ] = obj.GetMeshAreas(DT_tau);                                               % triangle areas
            [ DT_tau_triangle_edge_lengths ] = obj.GetMeshEdgeLengths(DT_tau);                              % edge lengths
            [ DT_tau_total_surface_integral ] = obj.GetTotalDensitySurfaceIntegral(DT_tau, density_tau);    % validate conservation

            obj.quality_indicators_history(end+1, :) = [min(DT_tau_trig_areas), max(DT_tau_trig_areas), min(DT_tau_triangle_edge_lengths), max(DT_tau_triangle_edge_lengths), DT_tau_total_surface_integral];
        end
        
        function VisualizeStep(obj, tstep)
            if(tstep < 0) %visualize last time step if not given or negative values
                tstep = obj.time_step_number;
            end
            
            obj.visualize_triangulation( obj.DT_history{tstep}, obj.density_history{tstep} );
            title(['Time - ', num2str(obj.time_history{tstep})]);
            xlim(obj.mesh_domain_limits(1, :)*2)
            ylim(obj.mesh_domain_limits(2, :))
            %caxis([0, max(max(cell2mat(density_history)))])
            caxis([0, 10]);
            colorpalette = colormap('jet');
            set(gca,'Color', colorpalette(1, :));   %assuming xaxis lower limit is zero, and density of areas outside of meshed surface is also zero
        end
        
        function VisualizeAnimation(obj, userdef_coupler, num_repeats)
            while(num_repeats ~= 0) %repeat num_repeats time, or non stop if negative
                
                for n = 1:obj.time_step_number
                    obj.VisualizeStep(n);
                    userdef_coupler.VisualizeStep(n);
                    
                    if(n < obj.time_step_number)
                        pause(obj.time_history{n+1} - obj.time_history{n})
                    end
                end
                
                num_repeats = num_repeats - 1;
            end
        end
        
        function SaveAnimation(obj, userdef_coupler, fname)
            h = figure;
            axis tight manual % this ensures that getframe() returns a consistent size
            filename = fname; %'../plots/may_31/sample_sim.gif';
            for n = 1:obj.time_step_number
                % Draw plot for y = x.^n

                obj.VisualizeStep(n);
                userdef_coupler.VisualizeStep(n);

                drawnow 
                  % Capture the plot as an image 
                  frame = getframe(h); 
                  im = frame2im(frame); 
                  [imind,cm] = rgb2ind(im,256); 
                  % Write to the GIF File 
                  if n == 1 
                      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
                  else 
                      imwrite(imind,cm,filename,'gif','WriteMode','append'); 
                  end 
              end
        end
    end
    
    methods(Static)
        function [ areas ] = GetMeshAreas(DT)
            trigs_x = reshape(DT.Points(DT.ConnectivityList', 1), [3, size(DT.ConnectivityList, 1)])';   %size = (num trigs, num vert per trig)
            trigs_y = reshape(DT.Points(DT.ConnectivityList', 2), [3, size(DT.ConnectivityList, 1)])';   %size = (num trigs, num vert per trig)
            areas = polyarea(trigs_x', trigs_y')';
        end
        
        function [ edge_lengths ] = GetMeshEdgeLengths( DT )
            all_edges = [DT.ConnectivityList(:, [1, 2]); DT.ConnectivityList(:, [2, 3]); DT.ConnectivityList(:, [1, 3])];
            unique_edges = unique(sort(all_edges, 2), 'rows');

            edge_end_1 = DT.Points(unique_edges(:, 1), :);
            edge_end_2 = DT.Points(unique_edges(:, 2), :);

            edge_lengths = sqrt(sum((edge_end_2 - edge_end_1).^2, 2));
        end
        
        function [ total_surface_integral ] = GetTotalDensitySurfaceIntegral( DT, density )
            total_surface_integral = density' * DiscreteTimeEvolvingMesh.GetMeshAreas( DT );
        end
        
        function [ centroids ] = GetMeshCentroids( DT )
            trigs_x = reshape(DT.Points(DT.ConnectivityList', 1), [3, size(DT.ConnectivityList, 1)])';   %size = (num trigs, num vert per trig)
            trigs_y = reshape(DT.Points(DT.ConnectivityList', 2), [3, size(DT.ConnectivityList, 1)])';   %size = (num trigs, num vert per trig)
            centroids = [mean(trigs_x, 2), mean(trigs_y, 2)];
        end
        
        function visualize_triangulation( DT_plot, density_plot )
            %ref: https://www.mathworks.com/matlabcentral/answers/165624-how-to-color-trisurf-faces

            trianglation_edge_alpha = 0;

            hh = trisurf(DT_plot.ConnectivityList, DT_plot.Points(:, 1), DT_plot.Points(:, 2), DT_plot.Points(:, 2)*0+1);
            set(gca,'CLim',[min(density_plot), max(density_plot)]);
            set(hh,'FaceColor','flat', 'FaceVertexCData',density_plot, 'CDataMapping','scaled');
            set(hh,'edgecolor', ([1 1 1])*0.3, 'EdgeAlpha', trianglation_edge_alpha)
            colorbar

            view(0,90)
        end
    end
    
end
