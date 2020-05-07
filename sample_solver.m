close all
clear 
clc


x_domain_lim = [0, 60];
v_domain_lim = [0, 20];


init_domain_p = [x_domain_lim(1), v_domain_lim(1); x_domain_lim(2), v_domain_lim(1); x_domain_lim(2), v_domain_lim(2); x_domain_lim(1), v_domain_lim(2); x_domain_lim(1), v_domain_lim(1)];   % (x, v)

[p,t,e] = pmesh(init_domain_p, 50, 6);
e=boundary_nodes(t);

v = zeros([size(p, 1), 1]);
v((p(:, 1) > 10 & p(:, 1) < 20 & p(:, 2) > 5 & p(:, 2) < 8)) = 1;

trisurf(t,p(:, 1),p(:, 2),v)

time = 0;       %time: seconds
dt = 0.1;       %t step size: seconds
T = 10;         %end time: seconds

p_t = p;
p_history = p_t;
t_index = 1;
for time = 0:dt:T
    p_tau = p_t;
    parfor p_ind = 1:size(p, 1)
        p_tau(p_ind, :) = [p_t(p_ind, 1) + p_t(p_ind, 2) * dt, p_t(p_ind, 2) + H(p_t, t, e, v, p_ind) * dt];
    end
    
    t_index = t_index + 1;
    p_t = p_tau;
    p_history(:, :, t_index) = p_t;
end




dx = 0.1;     %meter
dv = 0.1;     %meter/sec

x_list = x_domain_lim(1):dx:x_domain_lim(2);
v_list = v_domain_lim(1):dv:v_domain_lim(2);

[x_grid, y_grid] = meshgrid(x_list, v_list);
xy_domain = [x_grid(:), y_grid(:)];

vi = knnsearch(p_t, xy_domain);
Vq = v(vi);

figure
plot(sum(reshape(Vq, size(x_grid))))
%surf(x_grid, y_grid, reshape(Vq, size(x_grid)))


i = 1;
while(true)
    pause(0.1)
    vi = knnsearch(p_history(:, :, i), xy_domain);
    Vq = v(vi);
    %plot(sum(reshape(Vq, size(x_grid))))
    %surf(x_grid, y_grid, reshape(Vq, size(x_grid)))
    %draw on
    imagesc(reshape(Vq, size(x_grid)))
    set(gca,'YDir','normal')
    i = mod(i+1, size(p_history, 3)) + 1;
end