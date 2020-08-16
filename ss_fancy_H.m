function [ H_val ] = ss_fancy_H(density_mtx, boundaries_v_mtx, boundaries_t_mtx, eval_points)
%function [ H_val ] = H_IDM( evolving_mesh_obj, user_defined_coupler_obj, H_args_xv )
% DT ~ density surface mesh (currently implemented as delaunay triangulation)
% density ~ density value at inside each mesh element (triangle)
% H_args_xv ~ points at which to evaluate H function ~ each row is a 2d vector [x, v]
% assume_fixed_DT ~= 0 ==> triangulation is assumed not changing between
%                           iterations ==> this permits pre-computation of
%                           integral domains and relevant matrices
%                      ==> NOT IMPLEMENTED
    
    tmp_mtx = boundaries_v_mtx;
    flat_v_element_boundaries = [reshape(tmp_mtx(:, 1:end-1), size(tmp_mtx(:, 1:end-1), 1)*size(tmp_mtx(:, 1:end-1), 2), 1), reshape(tmp_mtx(:, 2:end), size(tmp_mtx(:, 2:end), 1)*size(tmp_mtx(:, 2:end), 2), 1)];

    tmp_mtx = boundaries_t_mtx;
    flat_t_element_boundaries = [reshape(tmp_mtx(:, 1:end-1), size(tmp_mtx(:, 1:end-1), 1)*size(tmp_mtx(:, 1:end-1), 2), 1), reshape(tmp_mtx(:, 2:end), size(tmp_mtx(:, 2:end), 1)*size(tmp_mtx(:, 2:end), 2), 1)];
    
    flat_density = reshape(density_mtx, size(density_mtx, 1) * size(density_mtx, 2), 1);
    density = flat_density;
    
    H_args_xv = eval_points;
    
    %----------------------------------
    flat_v_element_boundaries(isnan(flat_v_element_boundaries)) = 0;
    flat_t_element_boundaries(isnan(flat_t_element_boundaries)) = 0;
    density(isnan(density)) = 0;
    
    %DT = evolving_mesh_obj.DT;
    %density = evolving_mesh_obj.density;

    %------------- | parameters
    epsilon_nut = 1;
    d_nut = 2.5;
    V_max = 30;         %same as v_o for IDM??!
    alpha = 1;
    beta = 20;
    
    %H_IDM parameters       
    v_o = V_max;        %the velocity the vehicle would drive at in free traffic
    little_delta = 4;   
    s_o = 2;            %minimum spacings0:  a minimum desired net distance.  A car can not move if the distance from thecar in front is not at leasts0
    T = 1;              %desired time headwayT:  the minimum possible time to the vehicle in front (headway is the distancebetween vehicles in a transit system measured in time or space; the minimum headway is the shortestsuch distance or time achieved by a system without a reduction in the speed of vehicles)
    a = 1.3;              %accelerationa:  the maximum vehicle acceleration
    b = 2;              %comfortable braking decelerationb:  a positive number
    l = 4;              %gives the length of the vehicle
    
    
    %------------- | get mesh surface geometry information
    %NOTE!!!!!! area computation is approximate.. and assumes constant dt!!!! and assumes dt positive?! matters?
    [DT_centroids] = [mean(flat_v_element_boundaries, 2), mean(flat_t_element_boundaries, 2)]; %evolving_mesh_obj.GetMeshCentroids(DT);  %sorted by triangle id --> can pre-compute if DT is fixed
    [DT_areas] = diff(flat_v_element_boundaries')' * abs(nanmean(nanmean(diff(boundaries_t_mtx)))); %evolving_mesh_obj.GetMeshAreas(DT);          %sorted by triangle id --> can pre-compute if DT is fixed
    
	%------------- | process and evaluate IDM effects
    num_H_evaluation_points_t = size(H_args_xv, 1);
    num_DT_trigs_t = size(DT_centroids, 1);

    areas_t_mtx = repelem(DT_areas, 1, num_H_evaluation_points_t);                          %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)

    centroids_x_t_mtx = repelem(DT_centroids(:, 1), 1, num_H_evaluation_points_t);          %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)
    centroids_v_t_mtx = repelem(DT_centroids(:, 2), 1, num_H_evaluation_points_t);          %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)

    DT_x_t_mtx = repelem(H_args_xv(:, 1)', num_DT_trigs_t, 1);                 %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)
    DT_v_t_mtx = repelem(H_args_xv(:, 2)', num_DT_trigs_t, 1);                 %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)


    %H1_integral_domain = ( ((DT_x_t_mtx < centroids_x_t_mtx) & (centroids_x_t_mtx < DT_x_t_mtx+epsilon_nut)) ) & ( ((0 < centroids_v_t_mtx) & (centroids_v_t_mtx < max(DT.Points(:, 2)))) );


    y_mtx = centroids_x_t_mtx; %(H1_integral_domain);       %note: already incorporated into H1_integral_domain
    v_telda_mtx = centroids_v_t_mtx; %(H1_integral_domain); %note: already incorporated into H1_integral_domain
    integral_area_mtx = areas_t_mtx; %(H1_integral_domain);

    pre_computed_h_mtx = h(DT_x_t_mtx - y_mtx, epsilon_nut);
    %-------------
    
    %-------------
    H_IDM_evalpts_density_mtx = repelem(density, 1, num_H_evaluation_points_t);       %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)
    
	H_IDM_effect = a - a.* (H_args_xv(:, 2)./v_o).^little_delta - a.* sum(pre_computed_h_mtx .* H_IDM_evalpts_density_mtx .* ( (s_o + DT_v_t_mtx .* T + DT_v_t_mtx .* (DT_v_t_mtx - v_telda_mtx) ./ sqrt(2.*a.*b) ) ./ (y_mtx - DT_x_t_mtx - l) ).^2 .* integral_area_mtx)';
    
    H_IDM_effect(isnan(H_IDM_effect)) = 0;  %WARNING!! hard coded fix to numerical explosions caused by div by zero.. it should be ok for now, BUT fix later
    %-------------
    

    %------------- | evaluate total effect
    H_val = H_IDM_effect;
    %-------------



function [h_val] = h(x, epsilon_nut)
    h_val = exp(-1./( (epsilon_nut/2)^2 - (-x-epsilon_nut/2).^2 )) .* ((-epsilon_nut < x') & (x' < 0))';
    h_val(isnan(h_val)) = 0;


function [V_val] = V(x, d_nut, V_max)
    V_val = V_max * ((tanh(x/d_nut - 2) + tanh(2)) / (1+tanh(2)));


