function [ H_val ] = H( DT, DT_centroids, DT_areas, density, p_ind )
%H Summary of this function goes here
%   Detailed explanation goes here

    x_t = DT.Points(p_ind, 1);
    v_t = DT.Points(p_ind, 2);

    dx = 1;     %meter
    dv = 1;   %meter/sec
    
    epsilon_nut = 2;
    d_nut = 2;
    V_max = 20;
    
    %H1_sum_y = x_t:dx:x_t+epsilon_nut;
    %H1_sum_v_telda = 0:dv:max(DT.Points(:, 2));
    
    %[H1_sum_y_grid, H1_sum_v_telda_grid] = meshgrid(H1_sum_y, H1_sum_v_telda);
    %H1_integral_domain = [H1_sum_y_grid(:), H1_sum_v_telda_grid(:)];
    
    %trig_ids = DT.pointLocation(H1_integral_domain);
    %density_grid = zeros(size(trig_ids));
    %density_grid(~isnan(trig_ids)) = density(trig_ids(~isnan(trig_ids)));
    
    
    %[DT_centroids, DT_areas] = GetDelaunayCentroids(DT);  %sorted by triangle id
    H1_integral_domain = zeros([size(DT_centroids, 1), 1]);
    H1_integral_domain( ~((x_t < DT_centroids(:, 1)) & (DT_centroids(:, 1) < x_t+epsilon_nut)) ) = nan;
    H1_integral_domain( ~((0 < DT_centroids(:, 2)) & (DT_centroids(:, 2) < max(DT.Points(:, 2)))) ) = nan;
    
    density_grid = density;
    density_grid(isnan(H1_integral_domain)) = nan;
    
%     H1 = 0;
%     for y = H1_sum_y
%         for v_telda = H1_sum_v_telda
%             vi = knnsearch(p_t, [y, v_telda]);
%             Vq = v(vi);
%             
%             %[ti,bc] = pointLocation(t, [y, v_telda]);
%             %triVals = V(t(ti,:));
%             %Vq = dot(bc',triVals')';
%             
%             H1 = H1 + Vq * h(x_t - y, epsilon_nut) * (V(y - x_t, d_nut, V_max) - v_t);
%         end
%     end

    H1 = zeros([size(H1_integral_domain, 1), 1]);
    for integral_point_index = find(~isnan(H1_integral_domain))' %1:size(H1_integral_domain, 1)
        y = DT_centroids(integral_point_index, 1);
        v_telda = DT_centroids(integral_point_index, 2);
        integral_area = DT_areas(integral_point_index);
        H1(integral_point_index) = density_grid(integral_point_index) * h(x_t - y, epsilon_nut) * (V(y - x_t, d_nut, V_max) - v_t) * integral_area;
    end
    
    H_val = sum(H1);



function [h_val] = h(x, epsilon_nut)
    h_val = exp(-1/( (epsilon_nut/2)^2 - (-x-epsilon_nut/2)^2 )) * ((-epsilon_nut < x) & (x < 0));
    h_val(isnan(h_val)) = 0;


function [V_val] = V(x, d_nut, V_max)
    V_val = V_max * ((tanh(x/d_nut - 2) + tanh(2)) / (1+tanh(2)));


