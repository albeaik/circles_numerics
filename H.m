function [ H_val ] = H( DT, DT_centroids, DT_areas, density, assume_fixed_DT )
% assume_fixed_DT ~= 0 ==> triangulation is assumed not changing between
%                           iterations ==> this permits pre-computation of
%                           integral domains and relevant matrices
    
% be careful with global variables!!! especially pre_computation_done

    global pre_computation_done H1_integral_domain areas_t_mtx centroids_x_t_mtx centroids_v_t_mtx DT_x_t_mtx DT_v_t_mtx num_DT_points_t num_DT_trigs_t
    
    epsilon_nut = 2;
    d_nut = 2;
    V_max = 20;
    
    if(assume_fixed_DT == 0 | isempty(pre_computation_done)) %cold start
        num_DT_points_t = size(DT.Points, 1);
        num_DT_trigs_t = size(DT_centroids, 1);
        
        areas_t_mtx = repelem(DT_areas, 1, num_DT_points_t);         %size = (num trigs, num characteristic samples)

        centroids_x_t_mtx = repelem(DT_centroids(:, 1), 1, num_DT_points_t);         %size = (num trigs, num characteristic samples)
        centroids_v_t_mtx = repelem(DT_centroids(:, 2), 1, num_DT_points_t);         %size = (num trigs, num characteristic samples)

        DT_x_t_mtx = repelem(DT.Points(:, 1)', num_DT_trigs_t, 1);                 %size = (num trigs, num characteristic samples)
        DT_v_t_mtx = repelem(DT.Points(:, 2)', num_DT_trigs_t, 1);                 %size = (num trigs, num characteristic samples)

        H1_integral_domain = ( ((DT_x_t_mtx < centroids_x_t_mtx) & (centroids_x_t_mtx < DT_x_t_mtx+epsilon_nut)) ) & ( ((0 < centroids_v_t_mtx) & (centroids_v_t_mtx < max(DT.Points(:, 2)))) );
        
        pre_computation_done = 1;
    end
        
    density_t_mtx = repelem(density, 1, num_DT_points_t);         %size = (num trigs, num characteristic samples)
        
    y_mtx = centroids_x_t_mtx; %(H1_integral_domain);       %note: already incorporated into H1_integral_domain
    v_telda_mtx = centroids_v_t_mtx; %(H1_integral_domain); %note: already incorporated into H1_integral_domain
    integral_area_mtx = areas_t_mtx; %(H1_integral_domain);
    integral_density_mtx = density_t_mtx; %(H1_integral_domain);
    
    H_val = sum(H1_integral_domain .* integral_density_mtx .* h(DT_x_t_mtx - y_mtx, epsilon_nut) .* (V(y_mtx - DT_x_t_mtx, d_nut, V_max) - DT_v_t_mtx) .* integral_area_mtx);




function [h_val] = h(x, epsilon_nut)
    h_val = exp(-1./( (epsilon_nut/2)^2 - (-x-epsilon_nut/2).^2 )) .* ((-epsilon_nut < x') & (x' < 0))';
    h_val(isnan(h_val)) = 0;


function [V_val] = V(x, d_nut, V_max)
    V_val = V_max * ((tanh(x/d_nut - 2) + tanh(2)) / (1+tanh(2)));


