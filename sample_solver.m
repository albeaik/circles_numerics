close all
clear 
clc


x_domain_lim = [0, 60];
v_domain_lim = [0, 20];
n = 500;

init_domain_points = [x_domain_lim(1), v_domain_lim(1); x_domain_lim(2), v_domain_lim(1); x_domain_lim(2), v_domain_lim(2); x_domain_lim(1), v_domain_lim(2); x_domain_lim(1), v_domain_lim(1)];   % (x, v)
%init_domain_points = [init_domain_points; rand(n, 1)*x_domain_lim(2), rand(n, 1)*v_domain_lim(2)];

[p,t,e] = pmesh(init_domain_points, 50, 6);
%e=boundary_nodes(t);

DT.Points = p;
DT.ConnectivityList = t;

%DT = delaunayTriangulation(init_domain_points);
%p = DT.Points;
%t = DT.ConnectivityList;

%tmp_trig_centers = DT.incenter;
tmp_trig_centers = GetDelaunayCentroids( DT );
density = zeros([size(DT.ConnectivityList, 1), 1]);
density((tmp_trig_centers(:, 1) > 10 & tmp_trig_centers(:, 1) < 30 & tmp_trig_centers(:, 2) > 5 & tmp_trig_centers(:, 2) < 10)) = 1;

%tmp_trig_centers = DT.incenter;
tmp_trig_centers = GetDelaunayCentroids( DT );
%triplot(DT)
%trisurf(DT.ConnectivityList, DT.Points(:, 1), DT.Points(:, 2), density)
plot3(tmp_trig_centers(:, 1), tmp_trig_centers(:, 2), density, '.', 'markersize', 15)

time = 0;       %time: seconds
dt = 0.01;       %t step size: seconds
T = 10;          %end time: seconds

DT_t = DT;
DT_history{1} = DT_t;
density_t = density;
density_history{1} = density_t;
t_index = 1;
for time = 0:dt:T
    time
    DT_tau = DT_t;
    density_tau = density_t;
    
    [DT_t_centroids, DT_t_areas] = GetDelaunayCentroids(DT_t);  %sorted by triangle id
    
    DT_tau_points = DT_tau.Points;
    parfor point_ind = 1:size(DT.Points, 1)
        DT_tau_points(point_ind, :) = [DT_t.Points(point_ind, 1) + DT_t.Points(point_ind, 2) * dt, ...
                                    DT_t.Points(point_ind, 2) + H(DT_t, DT_t_centroids, DT_t_areas, density, point_ind) * dt];
    end
    DT_tau.Points = DT_tau_points;
    
    parfor trig_ind = 1:size(DT.ConnectivityList, 1)
        trig_t = DT_t.ConnectivityList(trig_ind, :);
        trig_tau = DT_tau.ConnectivityList(trig_ind, :);
        
        trig_t_coords = DT_t.Points(trig_t, :);
        trig_tau_coords = DT_tau.Points(trig_tau, :);
        
        density_tau(trig_ind) = (polyarea(trig_t_coords(:, 1), trig_t_coords(:, 2)) / polyarea(trig_tau_coords(:, 1), trig_tau_coords(:, 2))) * density_t(trig_ind);
    end
    
    t_index = t_index + 1;
    DT_t = DT_tau;
    DT_history{t_index} = DT_t;
    density_t = density_tau;
    density_history{t_index} = density_t;
end




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
    pause(0.01)
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
    
    i = mod(i+1, size(DT_history, 2)) + 1;
end