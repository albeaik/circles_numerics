close all
clear all
clc

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
T = 2.5;          %end time: seconds
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
    DT_tau_Points = [DT_t.Points(:, 1) + DT_t.Points(:, 2) * dt, ...
                                DT_t.Points(:, 2) + H(DT_t, density_t, Autonomous_state_t, DT_t.Points, assume_fixed_DT) * dt];
    
	% simulate the autonomous car ODE
    Autonomous_state_tau = Autonomous_state_t + [Autonomous_state_t(:, 2)*dt, (H(DT_t, density_t, Autonomous_state_t, Autonomous_state_t, assume_fixed_DT) + Autonomous_control_u) * dt];
                            
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
    xlim(x_domain_lim*2)
    ylim(v_domain_lim)
    %caxis([0, max(max(cell2mat(density_history)))])
    caxis([0, 10]);
    colorpalette = colormap('jet');
    set(gca,'Color', colorpalette(1, :));   %assuming xaxis lower limit is zero, and density of areas outside of meshed surface is also zero
    
    hold on
    plot3(Autonomous_state_history{end}(:, 1), Autonomous_state_history{end}(:, 2), 10*ones(size(Autonomous_state_history{end}(:, 2))), 'r*')
    hold off
    
    drawnow
end




% raw visualizations ---------------------------------------------------------


i = 1;
while(true)
    pause(0.1)
    i

    visualize_trig_trig( DT_history{i}, density_history{i} );
    title(['time - ', num2str(i/100)]);
    xlim(x_domain_lim)
    ylim(v_domain_lim)
    %caxis([0, max(max(cell2mat(density_history)))])
    caxis([0, 10])
    
    hold on
    plot3(Autonomous_state_history{i}(:, 1), Autonomous_state_history{i}(:, 2), 10*ones(size(Autonomous_state_history{i}(:, 2))), 'r*')
    hold off
    
    i = mod(i+1, size(DT_history, 2) -1 ) + 1;
end







% visualizations ---------------------------------------------------------
dx = 0.1;     %meter
dv = 0.1;     %meter/sec

x_list = x_domain_lim(1):dx:x_domain_lim(2);
v_list = v_domain_lim(1):dv:v_domain_lim(2);

[x_grid, y_grid] = meshgrid(x_list, v_list);
xy_domain = [x_grid(:), y_grid(:)];

% trig_ids = DT_t.pointLocation(xy_domain);
% density_grid = zeros(size(trig_ids));
% density_grid(~isnan(trig_ids)) = density_t(trig_ids(~isnan(trig_ids)));
% 
% figure
% plot(sum(reshape(density_grid, size(x_grid))))
%surf(x_grid, y_grid, reshape(Vq, size(x_grid)))


i = 1;
while(true)
    pause(0.1)
    i
    %vi = knnsearch(p_history(:, :, i), xy_domain);
    %Vq = v(vi);    
    
    %trig_ids = DT_history{i}.pointLocation(xy_domain);
    %density_grid = zeros(size(trig_ids));
    %density_grid(~isnan(trig_ids)) = density_history{i}(trig_ids(~isnan(trig_ids)));
    
    
    [DT_centroids, DT_areas] = GetDelaunayCentroids(DT_history{i});  %sorted by triangle id
    density_grid = griddata(DT_centroids(:, 1), DT_centroids(:, 2), density_history{i}, xy_domain(:, 1), xy_domain(:, 2));
    
    %plot(sum(reshape(Vq, size(x_grid))))
    %surf(x_grid, y_grid, reshape(Vq, size(x_grid)))
    %draw on
    
    imagesc(reshape(density_grid, size(x_grid)))
    colorbar
    caxis([0 10])
    set(gca,'YDir','normal')
    title(['time - ', num2str(i/100)])
    
    i = mod(i+10, size(DT_history, 2));
end





% visualizations ---------------------------------------------------------
figure
DT_plot = DT_tau;
triplot(DT_plot.ConnectivityList, DT_plot.Points(:, 1), DT_plot.Points(:, 2))

hold on
DT_plot = DT_tau_remished;
triplot(DT_plot.ConnectivityList, DT_plot.Points(:, 1), DT_plot.Points(:, 2), 'r')


%j = 55; figure; triplot(DT_history{j}.ConnectivityList, DT_history{j}.Points(:, 1), DT_history{j}.Points(:, 2))



% visualizations ---------------------------------------------------------
dx = 0.1;     %meter
dv = 0.1;     %meter/sec

x_list = x_domain_lim(1):dx:x_domain_lim(2);
v_list = v_domain_lim(1):dv:v_domain_lim(2);

[x_grid, y_grid] = meshgrid(x_list, v_list);
xy_domain = [x_grid(:), y_grid(:)];


DT_plot = DT_t;
density_plot = density_t;
[DT_centroids, DT_areas] = GetDelaunayCentroids(DT_plot);  %sorted by triangle id
density_grid = griddata(DT_centroids(:, 1), DT_centroids(:, 2), density_plot, xy_domain(:, 1), xy_domain(:, 2));

figure
imagesc(reshape(density_grid, size(x_grid)))
colorbar
caxis([0 10])
set(gca,'YDir','normal')


DT_plot = DT_tau;
density_plot = density_tau;
[DT_centroids, DT_areas] = GetDelaunayCentroids(DT_plot);  %sorted by triangle id
density_grid = griddata(DT_centroids(:, 1), DT_centroids(:, 2), density_plot, xy_domain(:, 1), xy_domain(:, 2));

figure
imagesc(reshape(density_grid, size(x_grid)))
colorbar
caxis([0 10])
set(gca,'YDir','normal')


DT_plot = DT_tau_remished;
density_plot = density_DT_tau_remished;
[DT_centroids, DT_areas] = GetDelaunayCentroids(DT_plot);  %sorted by triangle id
density_grid = griddata(DT_centroids(:, 1), DT_centroids(:, 2), density_plot, xy_domain(:, 1), xy_domain(:, 2));

figure
imagesc(reshape(density_grid, size(x_grid)))
colorbar
caxis([0 10])
set(gca,'YDir','normal')

%=====
DT_plot = DT_tau;
density_plot = density_tau;
[DT_centroids, DT_areas] = GetDelaunayCentroids(DT_plot);  %sorted by triangle id
density_grid_pre_remish = griddata(DT_centroids(:, 1), DT_centroids(:, 2), density_plot, xy_domain(:, 1), xy_domain(:, 2));

DT_plot = DT_tau_remished;
density_plot = density_DT_tau_remished;
[DT_centroids, DT_areas] = GetDelaunayCentroids(DT_plot);  %sorted by triangle id
density_grid_post_remish = griddata(DT_centroids(:, 1), DT_centroids(:, 2), density_plot, xy_domain(:, 1), xy_domain(:, 2));

figure
imagesc(reshape(density_grid_post_remish - density_grid_pre_remish, size(x_grid)))
colorbar
caxis([-0.1 0.1])
set(gca,'YDir','normal')



% visualizations ---------------------------------------------------------
DT_plot = DT_tau;
density_plot = density_tau;
[DT_centroids, DT_areas] = GetDelaunayCentroids(DT_plot);  %sorted by triangle id
density_grid = griddata(DT_centroids(:, 1), DT_centroids(:, 2), density_plot, xy_domain(:, 1), xy_domain(:, 2));

figure
imagesc(reshape(density_grid, size(x_grid)))
colorbar
caxis([0 10])
set(gca,'YDir','normal')

hold on
DT_plot = DT_tau;
triplot(DT_plot.ConnectivityList, DT_plot.Points(:, 1)*10, DT_plot.Points(:, 2)*10, 'g')

hold on
DT_plot = DT_tau_remished;
triplot(DT_plot.ConnectivityList, DT_plot.Points(:, 1)*10, DT_plot.Points(:, 2)*10, 'r')



%====
DT_plot = DT_tau_remished;
density_plot = density_DT_tau_remished;
[DT_centroids, DT_areas] = GetDelaunayCentroids(DT_plot);  %sorted by triangle id
density_grid = griddata(DT_centroids(:, 1), DT_centroids(:, 2), density_plot, xy_domain(:, 1), xy_domain(:, 2));

figure
imagesc(reshape(density_grid, size(x_grid)))
colorbar
caxis([0 10])
set(gca,'YDir','normal')

hold on
DT_plot = DT_tau;
triplot(DT_plot.ConnectivityList, DT_plot.Points(:, 1)*10, DT_plot.Points(:, 2)*10, 'g')

hold on
DT_plot = DT_tau_remished;
triplot(DT_plot.ConnectivityList, DT_plot.Points(:, 1)*10, DT_plot.Points(:, 2)*10, 'r')



%==============
%==============
return;
%============== | sample codes

%     %~~~~~~~~ | mesh/discretization maintanance
%     % (need to develop a strategy here)!! 1. evaluate, 2. adapt (in t, x, v), 3. propegate geneges or reset if necessary
%     % adaptation strategy!! 1. resampling, 2. retriangulation 3. interpolation
%     
%     %evaluate degeneracy (1. area based, 2. intersection based, .. etc)
%     refine_mesh_now = 0;
%     [DT_tau_trig_areas] = GetDelaunayAreas(DT_tau);
%     if((min(DT_tau_trig_areas) < mesh_subarea_minmax_limit(1)) || (max(DT_tau_trig_areas) > mesh_subarea_minmax_limit(2)))
%         dt = dt / 2;
%         if(dt < dt_minmax_limit(2))     %if time step is too small, reset dt and remesh instead
%             dt = dt_minmax_limit(1);
%             refine_mesh_now = 1;
%         else
%             continue;       %replace with a more robust reset strategy and ensure reset is proper
%         end
%     end
%------------    
%     %remesh   %constrained remesh and interpolation strategy.. good precision but code is junky and slow
%     if(refine_mesh_now)
%         DT_tau_remished = DT;                   %restart mesh to original.. and assume discontinuity constraints already incorporated into DT
%         DT_tau_remished.Points(unique(discontinuity_edge_list(:)), :) = DT_tau.Points(unique(discontinuity_edge_list(:)), :); %keep track of discontinuity location
%         DT_tau_remished_trig_centers = GetDelaunayCentroids( DT_tau_remished );
% 
%         density_DT_tau_remished = zeros([size(DT_tau_remished_trig_centers, 1), 1]);
%         %knn_point_id = knnsearch(DT_tau.Points, DT_tau_remished_trig_centers, 'K', 1);
% 
%         tmp_pt_set_1 = DT_tau.Points(DT_tau.ConnectivityList(:, 1), :);
%         tmp_pt_set_2 = DT_tau.Points(DT_tau.ConnectivityList(:, 2), :);
%         tmp_pt_set_3 = DT_tau.Points(DT_tau.ConnectivityList(:, 3), :);
% 
%         trig_edge_lengths = sqrt(sum([tmp_pt_set_1 - tmp_pt_set_2; tmp_pt_set_2 - tmp_pt_set_3; tmp_pt_set_1 - tmp_pt_set_2].^2, 2));
%         max_trig_edge_length = max(trig_edge_lengths);
% 
%         parfor i = 1:size(DT_tau_remished_trig_centers, 1)
%            %knn_pt_index = knn_point_id(i);
%            q_pt = DT_tau_remished_trig_centers(i, :);
% 
%            knn_pt_index = find(sqrt(sum((DT_tau.Points - repmat(q_pt, size(DT_tau.Points, 1), 1)).^2, 2)) <= 2*max_trig_edge_length);
% 
%            %tmp_inds = find(knn_pt_index == DT_tau.ConnectivityList(:));
%            %[tf, idx] = ismember(knn_pt_index, DT_tau.ConnectivityList(:)); %don't return repeat values
%            %tmp_inds = idx;
% 
%            tmp_inds = [];
%            for j = 1:size(knn_pt_index, 1)
%                tmp_inds = [tmp_inds; find(knn_pt_index(j) == DT_tau.ConnectivityList(:))];
%            end
% 
%            trig_ind = mod(tmp_inds-1, size(DT_tau.ConnectivityList, 1))+1;
% 
%            is_interior = zeros(size(trig_ind));
%            for j = 1:size(trig_ind, 1)
%                trig_xy = DT_tau.Points(DT_tau.ConnectivityList(trig_ind(j), :), :);
%                is_interior(j) = inpolygon(q_pt(1), q_pt(2), trig_xy(:, 1), trig_xy(:, 2));
%            end
% 
%            enclosing_trig = find(is_interior == 1); % could return multiple if point is at boundary
% 
%            if(~isempty(enclosing_trig))
%                density_DT_tau_remished(i) = density_tau(trig_ind(enclosing_trig(1))); %choose first if at boundary
%            else
%                density_DT_tau_remished(i) = 0;  %outside mesh.. extrapolate to zero
%            end
%         end
% 
%         %DT_tau_trig_centers = GetDelaunayCentroids( DT_tau );
%         %density_DT_tau_remished = griddata(DT_tau_trig_centers(:, 1), DT_tau_trig_centers(:, 2), density_tau, DT_tau_remished_trig_centers(:, 1), DT_tau_remished_trig_centers(:, 2));
% 
%         %if(sum(isnan(density_tau(:))) > 0)
%         %    disp('Warning: Nans from earlier steps are eliminated!!')
%         %end
% 
%         %density_DT_tau_remished(isnan(density_DT_tau_remished)) = 0; %Note: "griddata returns NaN for query points outside of the convex hull."
%         
%         assume_fixed_DT = 0;
%     else    %skip remishing
%         DT_tau_remished = DT_tau;
%         density_DT_tau_remished = density_tau;
%         assume_fixed_DT = 0;
%     end
    
%     % remesh back to initial mesh
%     if(refine_mesh_now)
%         DT_tau_remished = DT;       %new mish = origianl mish.. assuming DT wasn't changed since initialization
%         DT_tau_remished_trig_centers = GetDelaunayCentroids( DT_tau_remished );
%         DT_tau_trig_centers = GetDelaunayCentroids( DT_tau );
%         density_DT_tau_remished = griddata(DT_tau_trig_centers(:, 1), DT_tau_trig_centers(:, 2), density_tau, DT_tau_remished_trig_centers(:, 1), DT_tau_remished_trig_centers(:, 2));
%         
%         if(sum(isnan(density_tau(:))) > 0)
%             disp('Warning: Nans from earlier steps are eliminated!!')
%         end
%         
%         density_DT_tau_remished(isnan(density_DT_tau_remished)) = 0; %Note: "griddata returns NaN for query points outside of the convex hull."
%         assume_fixed_DT = 0;
%     else    %skip remishing
%         DT_tau_remished = DT_tau;
%         density_DT_tau_remished = density_tau;
%         assume_fixed_DT = 0;
%     end
%------------
%     % remesh: retriangulate DT_tau points
%     if(refine_mesh_now)
%         DT_tau_remished = DT;       %new mish = origianl mish.. assuming DT wasn't changed since initialization
%         DT_tau_remished.Points = DT_tau.Points; %constraints already set
%         DT_tau_remished_trig_centers = GetDelaunayCentroids( DT_tau_remished );
%         DT_tau_trig_centers = GetDelaunayCentroids( DT_tau );
%         density_DT_tau_remished = griddata(DT_tau_trig_centers(:, 1), DT_tau_trig_centers(:, 2), density_tau, DT_tau_remished_trig_centers(:, 1), DT_tau_remished_trig_centers(:, 2));
%         
%         if(sum(isnan(density_tau(:))) > 0)
%             disp('Warning: Nans from earlier steps are eliminated!!')
%         end
%         
%         density_DT_tau_remished(isnan(density_DT_tau_remished)) = 0; %Note: "griddata returns NaN for query points outside of the convex hull."
%         assume_fixed_DT = 0;
%     else    %skip remishing
%         DT_tau_remished = DT_tau;
%         density_DT_tau_remished = density_tau;
%         assume_fixed_DT = 0;
%     end