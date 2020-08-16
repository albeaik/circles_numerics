clear all
close all
clc


%---------------------| v discretization
v_init = 5;
v_end = 20;
dv = 0.1;      %used approximatly

num_v_bounderies = ceil((v_end - v_init) ./ dv);
v_init_boundaries = linspace(v_init, v_end, num_v_bounderies);

%---------------------| time discretization
T_init = 10;
T_end = 0;
dt_abs = 0.1;
dt_sign = sign(T_end-T_init);

num_t_snapshots = ceil(abs((T_end - T_init) / dt_abs));

%---------------------| pre-defined structures for solution surface
num_v_subdivs = num_v_bounderies - 1;
num_t_subdivs = num_t_snapshots;

%densities are assigned to subdivs
density_mtx = zeros([num_t_subdivs, num_v_subdivs])*nan;

%boundaries defining the subdivs
boundaries_v_mtx = zeros([num_t_subdivs, num_v_bounderies])*nan;
boundaries_t_mtx = zeros([num_t_subdivs, num_v_bounderies])*nan;

%subdiv v values
get_v_vals_mtx = @(mtx)  mean([mtx(1:end-1);mtx(2:end)]);

%get_v_vals_mtx_per_row = @(mtx) cellfun(get_v_vals_mtx, num2cell(mtx, 1)); %--> do tis

%---------------------| initial condition definitions
Q0_1 = @(t,v) (v >= 9) .* (v <= 11) + (v >= 12) .* (v <= 14) + (v >= 16) .* (v <= 17) * 0.5;
Q0_1 = @(t,v) sin(v).^2 .*v ;
Qt0_1 = @(t) t > 9;

Q0 = Q0_1;

t_init_condition = Qt0_1;

%---------------------
time = T_init;
dt = dt_sign * dt_abs;
solver_loop_index = 1;
while(dt_sign*(T_end - time) > 0)
    
    if(Qt0_1(time)) %use initial condition
        boundaries_v_mtx(solver_loop_index, :) = v_init_boundaries;
        boundaries_t_mtx(solver_loop_index, :) = zeros(size(v_init_boundaries)) + time;
        
        v_vals_mtx = get_v_vals_mtx(v_init_boundaries);
        density_mtx(solver_loop_index, :) = Q0(time, v_vals_mtx);
    else            %execute solver
        %reload
        old_boundaries_v_mtx = boundaries_v_mtx(solver_loop_index-1, :);
        old_boundaries_t_mtx = boundaries_t_mtx(solver_loop_index-1, :);
        old_density_mtx = density_mtx(solver_loop_index-1, :);
        
        %Compute! --> be careful when new v goes to zero given current sys dynamics
        new_boundaries_v_mtx = old_boundaries_v_mtx + ss_fancy_H(density_mtx./boundaries_v_mtx(:, 1:end-1), boundaries_v_mtx, boundaries_t_mtx, [(old_boundaries_t_mtx)', (old_boundaries_v_mtx)'])' ./ old_boundaries_v_mtx * dt; %dt should be signed? i.e., negative for backward soln
        new_boundaries_t_mtx = zeros(size(old_boundaries_t_mtx)) + time;
        
        new_density_vals = (diff(old_boundaries_v_mtx) ./ diff(new_boundaries_v_mtx)) .* old_density_mtx;
        
        %commit step
        boundaries_v_mtx(solver_loop_index, :) = new_boundaries_v_mtx;
        boundaries_t_mtx(solver_loop_index, :) = new_boundaries_t_mtx;
        density_mtx(solver_loop_index, :) = new_density_vals;
    end
    
    %-----| update solver time and index
    time = time + dt;
    solver_loop_index = solver_loop_index + 1;
end


%visualize
viz_v = zeros(size(density_mtx));
viz_t = zeros(size(density_mtx));
for i = 1:size(boundaries_v_mtx, 1)
    viz_v(i, :) = get_v_vals_mtx(boundaries_v_mtx(i, :));
    viz_t(i, :) = get_v_vals_mtx(boundaries_t_mtx(i, :));
end

figure
surf(viz_t, viz_v, density_mtx./viz_v)
colorbar
xlabel('Tau')
ylabel('v')
zlabel('density')