close all
clear all
clc

%load pre specified initial conditions
pre_specified_initial_conditions

%identify initial domain
x_domain_lim = [0, 20];
v_domain_lim = [0, 20];
x_resolution = 80;
v_resolution = 80;

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
T = 1;          %end time: seconds
remesh_cycle = 1;  %time steps
assume_fixed_DT = 0;

DT_t = DT;
DT_history{1} = DT_t;
density_t = density;
density_history{1} = density_t;
t_index = 1;
for time = 0:dt:T                   %NOTE: inside this loop, DT_.. changes from native native to non-native but compatible data type
    % prepare for this iteration
    time
    %DT_tau = DT_t;
    %density_tau = density_t;
    
    [DT_t_centroids] = GetDelaunayCentroids(DT_t);  %sorted by triangle id --> can pre-compute if DT is fixed
    [DT_t_areas] = GetDelaunayAreas(DT_t);  %sorted by triangle id --> can pre-compute if DT is fixed
    
    % simulate the characteristic equation
    DT_tau_Points = [DT_t.Points(:, 1) + DT_t.Points(:, 2) * dt, ...
                                DT_t.Points(:, 2) + H(DT_t, DT_t_centroids, DT_t_areas, density_t, assume_fixed_DT) * dt];
    
    %create DT_tau object ==> note matlab native DT object recomputes triangulation if points updated directly                       
	DT_tau.Points = DT_tau_Points;
    DT_tau.ConnectivityList = DT_t.ConnectivityList;
    
    % simulate density scaling
    [DT_t_trig_areas] = GetDelaunayAreas(DT_t);  %sorted by triangle id --> can pre-compute if DT is fixed
    [DT_tau_trig_areas] = GetDelaunayAreas(DT_tau);  %sorted by triangle id
    density_tau = (DT_t_trig_areas ./ DT_tau_trig_areas) .* density_t;
    
    %remesh
	DT_tau_remished = DT;                   %restart mesh to original.. and assume discontinuity constraints already incorporated into DT
    DT_tau_remished.Points(unique(discontinuity_edge_list(:)), :) = DT_tau.Points(unique(discontinuity_edge_list(:)), :); %keep track of discontinuity location
    DT_tau_remished_trig_centers = GetDelaunayCentroids( DT_tau_remished );
    DT_tau_trig_centers = GetDelaunayCentroids( DT_tau );
    density_DT_tau_remished = griddata(DT_tau_trig_centers(:, 1), DT_tau_trig_centers(:, 2), density_tau, DT_tau_remished_trig_centers(:, 1), DT_tau_remished_trig_centers(:, 2));

    if(sum(isnan(density_tau(:))) > 0)
        disp('Warning: Nans from earlier steps are eliminated!!')
    end

    density_DT_tau_remished(isnan(density_DT_tau_remished)) = 0; %Note: "griddata returns NaN for query points outside of the convex hull."

    
%     % remesh back to initial mesh
%     if(mod(time/dt, remesh_cycle) == 0)
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
%     else    %skip remishing
%         DT_tau_remished = DT_tau;
%         density_DT_tau_remished = density_tau;
%     end
    
    % prepare for next iteration
    t_index = t_index + 1;
    DT_t = DT_tau_remished; % DT_tau;
    DT_history{t_index} = DT_t;
    density_t = density_DT_tau_remished; % density_tau;
    density_history{t_index} = density_t;
end




% raw visualizations ---------------------------------------------------------


i = 1;
while(true)
    pause(0.1)
    i

    visualize_trig_trig( DT_history{i}, density_history{i} );
    title(['time - ', num2str(i/100)]);
    
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
