classdef DiscreteTimeEvolvingMesh < handle
    %DISCRETETIMEEVOLVINGMESH Summary of this class goes here
    %   Detailed explanation goes here
    %-- now assumes a triangulation.. can generalize to other types? 
    
    properties  %note: re-evaluate code... centralize all updates to these propertiese to ensure consistency
        DT;
        density;
        time;
        time_step_number;
        
        active_discontinuity_edges;     %--> note: these are implemented in ad hoc manner.. not consistent with history based commit structure.. can can cause issues if not carefully used
        direction_of_mesh_faces;        %--> note: these are implemented in ad hoc manner.. not consistent with history based commit structure.. can can cause issues if not carefully used
        
        mesh_domain_limits; %[x_domain_lim; v_domain_lim]   %2 by 2 matrix
        mesh_resolution;    %[x_resolution; v_resolution]   %2 by 1 vector
        
        time_history;
        DT_history;
        density_history;
        
        quality_indicators_history;
    end
    
    methods
        function obj = DiscreteTimeEvolvingMesh(initialization_option, initialization_parameters)
            %default constuctor function!
            %   initialization_option = {'oracle-function', 'pre-defined-initialization'}
            %   'oracle-function' -> initialization_parameters is a struct
            %                           initialization_parameters.Q0: queryable function such that init_density(x, v) <- Q0(x, v)
            %                           initialization_parameters.mesh_domain_limits [xmin, xmax; vmin, vmax]
            %                           initialization_parameters.mesh_resolution [number of mesh points x, number of mesh points v]
            %   'pre-defined-initialization' -> initialization_parameters is a struct
            %                           initialization_parameters.DT
            %                           initialization_parameters.density
            %old prototype: obj = DiscreteTimeEvolvingMesh(Q0, mesh_domain_limits, mesh_resolution)
            
            %set initial time
            time = 0;
            
            switch initialization_option
                case 'oracle-function'
                    [DT, density] = DiscreteTimeEvolvingMesh.GenerateMeshFromOracle(initialization_parameters);
                    mesh_domain_limits = initialization_parameters.mesh_domain_limits;
                    mesh_resolution = initialization_parameters.mesh_resolution;
                case 'pre-defined-initialization'
                    DT = initialization_parameters.DT;
                    density = initialization_parameters.density;
                    mesh_domain_limits = [min(DT.Points(:, 1)), max(DT.Points(:, 1)); min(DT.Points(:, 2)), max(DT.Points(:, 2))];
                    mesh_resolution = -1;
                otherwise
                    disp('Error!!')
            end

            %===========| construct object
            obj = obj@handle();
            
            obj.DT = DT;
            obj.density = density;
            obj.time = time;
            obj.time_step_number = 1;
            
            obj.active_discontinuity_edges = [];
            if(sum(strcmp(fieldnames(DT), 'Constraints')))  %consistent with matlab native object
                obj.active_discontinuity_edges = DT.Constraints;
            end
            
            obj.direction_of_mesh_faces = obj.GetTrigFaceDirections( DT );
            
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
        
        function CommitStep(obj, dt, newMesh, newMeshValues)
            obj.time_step_number = obj.time_step_number + 1;

            obj.time = obj.time + dt;
            obj.time_history{obj.time_step_number} = obj.time;

            if(sum(strcmp(fieldnames(newMesh), 'Constraints')))  %consistent with matlab native object
                obj.active_discontinuity_edges = newMesh.Constraints;           %only need to be updated when a remesh was conducted
            end
            obj.direction_of_mesh_faces = obj.GetTrigFaceDirections( newMesh ); %only need to be updated when a remesh was conducted

            obj.DT = newMesh;
            obj.DT_history{obj.time_step_number} = obj.DT;

            obj.density = newMeshValues;
            obj.density_history{obj.time_step_number} = obj.density;

            obj.MeshQualityIndicators(); %for now updates quality indicator array
        end
        
        function [newdt, newDT, newDensity, stepIsValid, deadEndFail] = DiscretizationValidationAndMaintanance(obj, dt, DT, density)
            %~~~~~~~~ | mesh/discretization maintanance
            % (need to develop a strategy here)!! 1. evaluate, 2. adapt (in t, x, v), 3. propegate geneges or reset if necessary
            % adaptation strategy!! 1. resampling, 2. retriangulation 3. interpolation
            % -----
            % temp implementation for now
            newdt = dt;
            newDT = DT;
            newDensity = density;
            stepIsValid = true;                     %new discretization is valid (either with or without maintanance)
            deadEndFail = false;                    %discretization is invaid and can't be fixed!
            
            min_dt_size = 1e-4;                 %HARD CODED VALUE!!!!!
            min_trig_area_threshold = 1e-4;     %HARD CODED VALUE!!!!!
            
            %-----------------------------------------------------------
            % assume face direction flip is a neccessary condition for triangule overlap 
            face_directions_now = obj.GetTrigFaceDirections( DT );
            face_directions_before = obj.direction_of_mesh_faces;
            triangleOverlap = (face_directions_now ~= face_directions_before);
            if(sum(triangleOverlap) > 0)
                stepIsValid = false; 
                newdt = dt / 2;
                
                
                if(newdt < min_dt_size)
                    deadEndFail = true;
                end
                
                return;         %NOTE: might need to revisit decision to return early from this function!!!!
            end
            
            %-----------------------------------------------------------
            [ DT_trig_areas ] = obj.GetMeshAreas(DT);                               % triangle areas   
            
            if(min(DT_trig_areas) < min_trig_area_threshold)                        %too small triangles
                disp('Remesh!!!')
                
                %-----------------| regenerate mesh with discontinuity constraints
                discont_edges = obj.active_discontinuity_edges;
                discont_points = [];
                if(~isempty(discont_edges)) %extract points given references from edges, and remap indecies to isolated point set
                    [unique_disc_pt_indecies, ia, ic] = unique([discont_edges(:, 1); discont_edges(:, 2)]);
                    disc_pts = DT.Points(unique_disc_pt_indecies, :);
                    new_pt_indecies = 1:size(unique_disc_pt_indecies, 1);
                    new_discont_edges_flat = new_pt_indecies(ic);
                    new_discont_edges = reshape(new_discont_edges_flat, [size(new_discont_edges_flat, 2)/2, 2]);
                    
                    discont_edges = new_discont_edges;
                    discont_points = disc_pts;
                end
                
                %update domain limits
                mesh_domain_limits = [min(DT.Points(:, 1)), max(DT.Points(:, 1)); min(DT.Points(:, 2)), max(DT.Points(:, 2))];
                
                %sample domain
                x=linspace(mesh_domain_limits(1, 1), mesh_domain_limits(1, 2), obj.mesh_resolution(1));
                v=linspace(mesh_domain_limits(2, 1), mesh_domain_limits(2, 2), obj.mesh_resolution(2));
                [x_grid, v_grid]=meshgrid(x, v);

                init_domain_points = [x_grid(:), v_grid(:)];
                init_domain_points = [discont_points; init_domain_points];

                if(~isempty(discont_edges))
                    newDT = delaunayTriangulation(init_domain_points, discont_edges);
                else
                    newDT = delaunayTriangulation(init_domain_points);
                end
                
                %-----------------| interpolate new density values
                new_trig_centers = DiscreteTimeEvolvingMesh.GetMeshCentroids( newDT );
                newDensity = zeros([size(newDT.ConnectivityList, 1), 1]);
                
                %need to improve computational efficiency!!
                parfor i = 1:size(newDT.ConnectivityList, 1)   %loop over all new triangles
                    for j = 1:size(DT.ConnectivityList, 1)  %loop over all old triangles
                        q_pt = new_trig_centers(i, :);
                        trig_xy = DT.Points(DT.ConnectivityList(j, :), :);
                        if(inpolygon(q_pt(1), q_pt(2), trig_xy(:, 1), trig_xy(:, 2)))
                            newDensity(i) = density(j);
                            continue;
                        end
                    end
                end
                
               %-----------------| plot
               figure
               subplot(1, 2, 1)
               triplot(DT.ConnectivityList, DT.Points(:, 1), DT.Points(:, 2))
               subplot(1, 2, 2)
               triplot(newDT.ConnectivityList, newDT.Points(:, 1), newDT.Points(:, 2))
            end
        end
        
        function MeshQualityIndicators(obj)
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
        
        function VisualizeMultipleSteps(obj, userdef_coupler, plt_times)
            %fighandle = figure
            %plt_times = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3];
            dt = obj.time_history{2}-obj.time_history{1};     %assume constant dt
            for t_ind = 1:size(plt_times, 2)
                subplot(1, size(plt_times, 2), t_ind)
                time = plt_times(t_ind);
                n = plt_times(t_ind)/dt+1; %find(cell2mat(solutions{index}.q.time_history) == time);
                obj.VisualizeStep(n);
                userdef_coupler.VisualizeStep(n);

                %xlim([0, 15])
                %ylim([5, 12])
                %caxis([0, 2]);
                %colorbar off;

                %if(t_ind == 1)
                %    ylabel('Characteristic Solver')
                %end
                %xlabel('t')
                %set(gca,'FontSize',15);
            end
            %set(gcf, 'Position',  [74         614        1562         334])
            set(gcf, 'Position',  [1         760        1680         188])
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
        function [DT, density] = GenerateMeshFromOracle(parameters)
            %parse
            Q0 = parameters.Q0;
            mesh_domain_limits = parameters.mesh_domain_limits;
            mesh_resolution = parameters.mesh_resolution;
            
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
            
            %process Q0 and get discontinuity
            discont_points = [];
            discont_edges = [];
            if(isa(Q0,'function_handle'))       %for backward compatibility
                Q0_oracle = Q0;
            else                                %implementation july 14, 2020
                Q0_oracle = Q0.oracle;
                if(~isempty(Q0.discontinuity_boundary_points) && ~isempty(Q0.discontinuity_boundary_edges))
                    discont_points = Q0.discontinuity_boundary_points;
                    discont_edges = Q0.discontinuity_boundary_edges;
                end
            end
            
            %sample domain
            x=linspace(mesh_domain_limits(1, 1), mesh_domain_limits(1, 2), mesh_resolution(1));
            v=linspace(mesh_domain_limits(2, 1), mesh_domain_limits(2, 2), mesh_resolution(2));
            [x_grid, v_grid]=meshgrid(x, v);
            
            init_domain_points = [x_grid(:), v_grid(:)];
            init_domain_points = [discont_points; init_domain_points];

            if(~isempty(discont_edges))
                DT = delaunayTriangulation(init_domain_points, discont_edges);
            else
                DT = delaunayTriangulation(init_domain_points);
            end
            
            %p = DT.Points;
            %t = DT.ConnectivityList;
                
            %--- generate density profile
            %tmp_trig_centers = DT.incenter;
            tmp_trig_centers = DiscreteTimeEvolvingMesh.GetMeshCentroids( DT );
            density = zeros([size(DT.ConnectivityList, 1), 1]);
            density(Q0_oracle(tmp_trig_centers(:, 1), tmp_trig_centers(:, 2)) == 1) = 1;

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
        end
        
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
        
        function [trig_face_direction] = GetTrigFaceDirections( DT )
            %comptes direction of faces assuming trianguler grid in 2D
            %using corss product approach.
            
            Pt_tail = DT.Points(DT.ConnectivityList(:, 1), :);
            Pt_head_1 = DT.Points(DT.ConnectivityList(:, 2), :);
            Pt_head_2 = DT.Points(DT.ConnectivityList(:, 3), :);
            
            trig_vect_1 = [Pt_head_1 - Pt_tail, zeros([size(Pt_head_1, 1), 1])];
            trig_vect_2 = [Pt_head_2 - Pt_tail, zeros([size(Pt_head_1, 1), 1])];
            
            normals_to_trigs = cross(trig_vect_1, trig_vect_2, 2);
            
            trig_face_direction = sign(normals_to_trigs(:, 3));
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

