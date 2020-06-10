function [ H_val ] = H( DT, density, H_args_xv, assume_fixed_DT )
% DT ~ density surface mesh (currently implemented as delaunay triangulation)
% density ~ density value at inside each mesh element (triangle)
% H_args_xv ~ points at which to evaluate H function ~ each row is a 2d vector [x, v]
% assume_fixed_DT ~= 0 ==> triangulation is assumed not changing between
%                           iterations ==> this permits pre-computation of
%                           integral domains and relevant matrices
    
% be careful with global variables!!! especially pre_computation_done
% be careful with precomputations!!! not fully implemented

    global pre_computation_done H1_integral_domain areas_t_mtx centroids_x_t_mtx centroids_v_t_mtx DT_x_t_mtx DT_v_t_mtx num_H_evaluation_points_t num_DT_trigs_t integral_area_mtx pre_computed_h_mtx pre_computed_Vterm_mtx
    global sparse_H1_integral_domain sparse_pre_computed_h_mtx sparse_pre_computed_Vterm_mtx sparse_integral_area_mtx
    global sparse_h_V_area_terms_H1 sparse_h_V_area_terms_H2
    
    epsilon_nut = 1;
    d_nut = 2.5;
    V_max = 20;
    beta = 1;
    
    [DT_centroids] = GetDelaunayCentroids(DT);  %sorted by triangle id --> can pre-compute if DT is fixed
    [DT_areas] = GetDelaunayAreas(DT);          %sorted by triangle id --> can pre-compute if DT is fixed
    
    if(assume_fixed_DT == 0 | isempty(pre_computation_done)) %cold start
        num_H_evaluation_points_t = size(H_args_xv, 1);
        num_DT_trigs_t = size(DT_centroids, 1);
        
        areas_t_mtx = repelem(DT_areas, 1, num_H_evaluation_points_t);                          %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)

        centroids_x_t_mtx = repelem(DT_centroids(:, 1), 1, num_H_evaluation_points_t);          %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)
        centroids_v_t_mtx = repelem(DT_centroids(:, 2), 1, num_H_evaluation_points_t);          %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)

        DT_x_t_mtx = repelem(H_args_xv(:, 1)', num_DT_trigs_t, 1);                 %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)
        DT_v_t_mtx = repelem(H_args_xv(:, 2)', num_DT_trigs_t, 1);                 %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)

        
        H1_integral_domain = ( ((DT_x_t_mtx < centroids_x_t_mtx) & (centroids_x_t_mtx < DT_x_t_mtx+epsilon_nut)) ) & ( ((0 < centroids_v_t_mtx) & (centroids_v_t_mtx < max(DT.Points(:, 2)))) );
        
        
        y_mtx = centroids_x_t_mtx; %(H1_integral_domain);       %note: already incorporated into H1_integral_domain
        v_telda_mtx = centroids_v_t_mtx; %(H1_integral_domain); %note: already incorporated into H1_integral_domain
        integral_area_mtx = areas_t_mtx; %(H1_integral_domain);

        pre_computed_h_mtx = h(DT_x_t_mtx - y_mtx, epsilon_nut);
        pre_computed_Vterm_mtx = (V(y_mtx - DT_x_t_mtx, d_nut, V_max) - DT_v_t_mtx);
        
        %sparcify
        %sparse_H1_integral_domain = sparse(H1_integral_domain);
        %sparse_pre_computed_h_mtx = sparse(sparse_H1_integral_domain .* pre_computed_h_mtx);
        %sparse_pre_computed_Vterm_mtx = sparse(sparse_H1_integral_domain .* pre_computed_Vterm_mtx);
        %sparse_integral_area_mtx = sparse(sparse_H1_integral_domain .* integral_area_mtx);
        
        sparse_h_V_area_terms_H1 = sparse(H1_integral_domain .* pre_computed_h_mtx .* pre_computed_Vterm_mtx .* integral_area_mtx);
        sparse_h_V_area_terms_H2 = sparse(pre_computed_h_mtx .* (DT_v_t_mtx - v_telda_mtx) .* integral_area_mtx);
        
        pre_computation_done = 1;
    end
    
    H1_val = sparse_h_V_area_terms_H1' * density;
    H2_val = - beta .* sparse_h_V_area_terms_H2' * density;
    H_val = H1_val + H2_val;
        
%     density_t_mtx = repelem(density, 1, num_DT_points_t);         %size = (num trigs, num characteristic samples)
%     integral_density_mtx = density_t_mtx; %(H1_integral_domain);
%     sparse_integral_density_mtx = sparse(sparse_H1_integral_domain .* integral_density_mtx);
    
%     H_val = sum(sparse_H1_integral_domain .* sparse_integral_density_mtx .* sparse_pre_computed_h_mtx .* sparse_pre_computed_Vterm_mtx .* sparse_integral_area_mtx);




function [h_val] = h(x, epsilon_nut)
    h_val = exp(-1./( (epsilon_nut/2)^2 - (-x-epsilon_nut/2).^2 )) .* ((-epsilon_nut < x') & (x' < 0))';
    h_val(isnan(h_val)) = 0;


function [V_val] = V(x, d_nut, V_max)
    V_val = V_max * ((tanh(x/d_nut - 2) + tanh(2)) / (1+tanh(2)));


