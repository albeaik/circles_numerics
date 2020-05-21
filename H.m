function [ H_val ] = H( DT, DT_centroids, DT_areas, density, p_ind )
%H Summary of this function goes here
%   Detailed explanation goes here

    x_t = DT.Points(p_ind, 1);
    v_t = DT.Points(p_ind, 2);

    dx = 1;     %meter
    dv = 1;     %meter/sec
    
    epsilon_nut = 2;
    d_nut = 2;
    V_max = 20;
    
    H1_integral_domain = ( ((x_t < DT_centroids(:, 1)) & (DT_centroids(:, 1) < x_t+epsilon_nut)) ) & ( ((0 < DT_centroids(:, 2)) & (DT_centroids(:, 2) < max(DT.Points(:, 2)))) );

    y = DT_centroids((H1_integral_domain), 1);
    v_telda = DT_centroids((H1_integral_domain), 2);
    integral_area = DT_areas((H1_integral_domain));
    H_val = sum(density((H1_integral_domain)) .* h(x_t - y, epsilon_nut)' .* (V(y - x_t, d_nut, V_max) - v_t) .* integral_area);




function [h_val] = h(x, epsilon_nut)
    h_val = exp(-1/( (epsilon_nut/2)^2 - (-x-epsilon_nut/2).^2 )) .* ((-epsilon_nut < x') & (x' < 0));
    h_val(isnan(h_val)) = 0;


function [V_val] = V(x, d_nut, V_max)
    V_val = V_max * ((tanh(x/d_nut - 2) + tanh(2)) / (1+tanh(2)));


